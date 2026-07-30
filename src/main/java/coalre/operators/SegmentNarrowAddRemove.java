package coalre.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Add-Remove operator that performs narrow exchange-like moves on individual segment trees.
 *
 * ADD: At a coalescent event for segment s, add a reassortment that diverts segment s
 *      from one child lineage to the uncle lineage, creating a narrow-exchange topology
 *      for that segment only.
 *
 * REMOVE: Remove a single-segment reassortment edge that was added by this operation,
 *         restoring the original coalescent structure for that segment.
 *
 * This maintains MCMC reversibility by ensuring the add and remove moves are exact inverses.
 */
public class SegmentNarrowAddRemove extends DivertSegmentOperator {

    public Input<Double> addProbabilityInput = new Input<>("addProbability",
            "Probability of adding a reassortment edge (vs removing)", 0.5);

    public Input<Double> sourceTimeScaleInput = new Input<>("sourceTimeScale",
            "Scale parameter for sampling source time relative to edge length", 0.5);

    private double addProb;
    private double sourceTimeScale;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        addProb = addProbabilityInput.get();
        sourceTimeScale = sourceTimeScaleInput.get();
    }

    @Override
    public double networkProposal() {
        double logHR = 0.0;

        if (Randomizer.nextDouble() < addProb) {
            logHR = addNarrowReassortment();
            logHR += Math.log((1.0 - addProb) / addProb);
        } else {
            logHR = removeNarrowReassortment();
            logHR += Math.log(addProb / (1.0 - addProb));
        }

        return logHR;
    }

    /**
     * ADD MOVE: Add a reassortment edge that creates a narrow-exchange topology
     * for a single segment tree.
     *
     * Works at the TREE level:
     * 1. Pick a segment tree
     * 2. Select grandparent node in tree (with internal children)
     * 3. Pick which child of grandparent to swap with uncle
     * 4. Add reassortment in network that swaps this child with uncle FOR THIS SEGMENT ONLY
     * 5. This creates narrow exchange topology in the segment tree
     */
    private double addNarrowReassortment() {
        double logHR = 0.0;

        // Step 1: Pick a random segment
        int segment = Randomizer.nextInt(network.getSegmentCount());

        // Step 2: Get a coalescent node in the segment tree (grandparent node)
        final int internalNodes = segmentTrees.get(segment).getInternalNodeCount();
        if (internalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        // Pick a random internal node to be the grandparent
        Node grandParent = segmentTrees.get(segment).getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));

        // Ensure grandparent has at least one internal child (for narrow exchange)
        while (grandParent.getLeft().isLeaf() && grandParent.getRight().isLeaf()) {
            grandParent = segmentTrees.get(segment).getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        }

        // Count valid grandparents for HR
        int validGP = 0;
        for (int i = internalNodes + 1; i < 1 + 2 * internalNodes; i++) {
            validGP += isg(segmentTrees.get(segment).getNode(i));
        }
        logHR -= Math.log(1.0 / validGP);

        // Step 3: Identify parent and uncle in the TREE
        Node parentNode = grandParent.getLeft();
        Node uncleNode = grandParent.getRight();

        // Parent should be higher (more recent) than uncle
        if (parentNode.getHeight() < uncleNode.getHeight()) {
            parentNode = grandParent.getRight();
            uncleNode = grandParent.getLeft();
        }

        if (parentNode.isLeaf()) {
            return Double.NEGATIVE_INFINITY;
        }

        // Step 4: Pick which child of parent to swap with uncle
        Node childToSwap = parentNode.getLeft();
        if (Randomizer.nextBoolean()) {
            childToSwap = parentNode.getRight();
        }
        logHR -= Math.log(0.5);

        // Now we need to perform the swap in the NETWORK for this segment only
        // The tree currently shows: grandParent -> parent -> (childToSwap, otherChild)
        //                                        -> uncle
        //
        // After swap, tree should show: grandParent -> parent -> (uncle, otherChild)
        //                                            -> childToSwap
        //
        // To do this in the network, we add reassortment that diverts segment from:
        // - childToSwap's lineage to uncle's lineage at grandParent height
        // This makes segment see childToSwap and uncle swapped

        // Step 5: Find network edges corresponding to these tree nodes
        int grandParentNr = grandParent.getNr();

        // Find network edge where grandParent coalesces
        List<NetworkEdge> gpEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.parentNode.isCoalescence())
                .filter(e -> e.parentNode.segmentIndices != null)
                .filter(e -> e.parentNode.segmentIndices[segment] == grandParentNr)
                .filter(e -> e.parentNode.getChildEdges().get(0).hasSegments.get(segment))
                .filter(e -> e.parentNode.getChildEdges().get(1).hasSegments.get(segment))
                .collect(Collectors.toList());

        if (gpEdges.size() != 2) {
            return Double.NEGATIVE_INFINITY;
        }

        NetworkNode gpNetworkNode = gpEdges.get(0).parentNode;
        double grandParentHeight = gpNetworkNode.getHeight();

        // Find which network edge corresponds to parent node in tree
        int parentNr = parentNode.getNr();
        NetworkEdge parentNetworkEdge = null;
        NetworkEdge uncleNetworkEdge = null;

        for (NetworkEdge e : gpNetworkNode.getChildEdges()) {
            if (e.hasSegments.get(segment)) {
                int nodeNr = getSegmentTreeNodeNr(e, segment);
                if (nodeNr == parentNr) {
                    parentNetworkEdge = e;
                } else {
                    uncleNetworkEdge = e;
                }
            }
        }

        if (parentNetworkEdge == null || uncleNetworkEdge == null) {
            return Double.NEGATIVE_INFINITY;
        }

        // Find network edge corresponding to childToSwap
        int childNr = childToSwap.getNr();
        NetworkEdge childNetworkEdge = findNetworkEdgeForTreeNode(parentNetworkEdge, childNr, segment);

        if (childNetworkEdge == null) {
            return Double.NEGATIVE_INFINITY;
        }

        // Step 6: Sample source time on child's network edge
        double minSourceTime = childNetworkEdge.childNode.getHeight();
        double maxSourceTime = Math.min(childNetworkEdge.parentNode.getHeight(), grandParentHeight);

        if (maxSourceTime <= minSourceTime) {
            return Double.NEGATIVE_INFINITY;
        }

        double sourceTime = minSourceTime + Randomizer.nextDouble() * (maxSourceTime - minSourceTime) * sourceTimeScale;
        logHR -= Math.log(sourceTimeScale / (maxSourceTime - minSourceTime));

        // Step 7: Find destination edge on uncle's lineage at grandParent height
        NetworkEdge destEdge = getEdgeWithHeight(uncleNetworkEdge, grandParentHeight, segment);

        if (destEdge == null) {
            return Double.NEGATIVE_INFINITY;
        }

        double destTime = grandParentHeight;

        // Step 8: Add the reassortment edge that diverts this segment
        logHR += addSegmentDivertingReassortment(childNetworkEdge, sourceTime,
                                                  destEdge, destTime, segment);

        if (logHR == Double.NEGATIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        }

        // Step 9: Calculate reverse probability
        int nRemovableEdges = countRemovableNarrowReassortments();

        if (nRemovableEdges == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        logHR += Math.log(1.0 / nRemovableEdges);

        return logHR;
    }

    /**
     * REMOVE MOVE: Remove a single-segment reassortment edge that creates
     * a narrow-exchange topology, restoring the coalescent structure.
     */
    private double removeNarrowReassortment() {
        double logHR = 0.0;

        // Step 1: Find all removable reassortment edges
        List<NetworkEdge> removableEdges = findRemovableNarrowReassortments();

        if (removableEdges.isEmpty()) {
            return Double.NEGATIVE_INFINITY;
        }

        // Step 2: Select an edge to remove
        NetworkEdge edgeToRemove = removableEdges.get(Randomizer.nextInt(removableEdges.size()));
        logHR -= Math.log(1.0 / removableEdges.size());

        // Step 3: Extract parameters before removal (for reverse probability)
        int segment = edgeToRemove.hasSegments.nextSetBit(0);
        NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
        NetworkNode grandParentNode = edgeToRemove.parentNode;

        // Before removal, get the tree node numbers involved
        int gpTreeNodeNr = (grandParentNode.segmentIndices != null) ?
                           grandParentNode.segmentIndices[segment] : -1;

        // Step 4: Remove the reassortment edge
        logHR += removeReassortmentEdge(edgeToRemove);

        if (logHR == Double.NEGATIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        }

        // Step 5: Calculate reverse probability (probability of adding this edge back)
        // This must match the ADD operation exactly!

        // After removal, the segment tree has been restored to coalescent structure
        // Count valid grandparents in the segment tree (same as ADD operation)
        final int internalNodes = segmentTrees.get(segment).getInternalNodeCount();

        if (internalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        int validGPafter = 0;
        for (int i = internalNodes + 1; i < 1 + 2 * internalNodes; i++) {
            validGPafter += isg(segmentTrees.get(segment).getNode(i));
        }

        // Find the grandparent node in the tree that corresponds to where we removed from
        Node grandParentTree = null;
        if (gpTreeNodeNr >= 0) {
            try {
                grandParentTree = segmentTrees.get(segment).getNode(gpTreeNodeNr);
            } catch (Exception e) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        if (grandParentTree == null || grandParentTree.isLeaf()) {
            return Double.NEGATIVE_INFINITY;
        }

        // Reverse probability for selecting this grandparent
        logHR += Math.log(1.0 / validGPafter);

        // Reverse probability for choosing which child to swap (was 0.5)
        logHR += Math.log(0.5);

        // Reverse probability for sampling source time
        // Need to reconstruct the edge that was the source
        NetworkEdge reconstructedSourceEdge = sourceEdge;
        double minSourceTime = reconstructedSourceEdge.childNode.getHeight();
        double maxSourceTime = Math.min(reconstructedSourceEdge.parentNode.getHeight(), grandParentNode.getHeight());

        if (maxSourceTime <= minSourceTime) {
            return Double.NEGATIVE_INFINITY;
        }

        logHR += Math.log(sourceTimeScale / (maxSourceTime - minSourceTime));

        return logHR;
    }

    /**
     * Add a reassortment edge that diverts a single segment from source to dest
     */
    private double addSegmentDivertingReassortment(NetworkEdge sourceEdge, double sourceTime,
                                                     NetworkEdge destEdge, double destTime,
                                                     int segment) {
        double logHR = 0.0;

        // Create source node
        NetworkNode sourceNode = new NetworkNode();
        sourceNode.setHeight(sourceTime);

        NetworkNode oldSourceParent = sourceEdge.parentNode;
        oldSourceParent.removeChildEdge(sourceEdge);
        sourceNode.addChildEdge(sourceEdge);

        NetworkEdge newSourceParentEdge = new NetworkEdge();
        sourceNode.addParentEdge(newSourceParentEdge);
        oldSourceParent.addChildEdge(newSourceParentEdge);
        newSourceParentEdge.hasSegments = (BitSet) sourceEdge.hasSegments.clone();

        // Create dest node
        NetworkNode destNode = new NetworkNode();
        destNode.setHeight(destTime);

        NetworkNode oldDestParent = destEdge.parentNode;
        oldDestParent.removeChildEdge(destEdge);
        destNode.addChildEdge(destEdge);

        NetworkEdge newDestParentEdge = new NetworkEdge();
        destNode.addParentEdge(newDestParentEdge);
        oldDestParent.addChildEdge(newDestParentEdge);
        newDestParentEdge.hasSegments = (BitSet) destEdge.hasSegments.clone();

        // Create reassortment edge
        NetworkEdge reassortEdge = new NetworkEdge();
        sourceNode.addParentEdge(reassortEdge);
        destNode.addChildEdge(reassortEdge);
        reassortEdge.hasSegments = new BitSet();

        // Divert only the selected segment
        BitSet segmentToDivert = new BitSet();
        segmentToDivert.set(segment);

        logHR += divertSegments(reassortEdge, newSourceParentEdge, segmentToDivert);

        // Update networkEdges list
        networkEdges.add(newSourceParentEdge);
        networkEdges.add(newDestParentEdge);
        networkEdges.add(reassortEdge);

        return logHR;
    }

    /**
     * Remove a reassortment edge (standard removal)
     */
    private double removeReassortmentEdge(NetworkEdge edgeToRemove) {
        double logHR = 0.0;

        NetworkNode nodeToRemove = edgeToRemove.childNode;
        NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
        NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;

        networkEdges.remove(edgeToRemove);
        networkEdges.remove(edgeToRemoveSpouse);

        // Divert segments back
        BitSet segsToDivert = (BitSet) edgeToRemove.hasSegments.clone();
        logHR += divertSegments(edgeToRemoveSpouse, edgeToRemove, segsToDivert);

        // For single segment, no subset probability needed (deterministic)
        // logHR += 0.0 for selecting 1 segment out of 1 segment

        // Remove nodes
        NetworkEdge edgeToExtend = nodeToRemove.getChildEdges().get(0);
        nodeToRemove.removeChildEdge(edgeToExtend);
        nodeToRemove.removeParentEdge(edgeToRemove);
        nodeToRemove.removeParentEdge(edgeToRemoveSpouse);
        edgeToRemoveSpouseParent.removeChildEdge(edgeToRemoveSpouse);
        edgeToRemoveSpouseParent.addChildEdge(edgeToExtend);

        NetworkNode secondNodeToRemove = edgeToRemove.parentNode;
        NetworkEdge secondEdgeToExtend = getSisterEdge(edgeToRemove);

        secondNodeToRemove.removeChildEdge(secondEdgeToExtend);
        secondNodeToRemove.removeChildEdge(edgeToRemove);

        networkEdges.remove(secondNodeToRemove.getParentEdges().get(0));

        if (secondNodeToRemove.getParentEdges().get(0).isRootEdge()) {
            network.setRootEdge(secondEdgeToExtend);
        } else {
            NetworkEdge secondNodeToRemoveParentEdge = secondNodeToRemove.getParentEdges().get(0);
            NetworkNode secondNodeToRemoveParent = secondNodeToRemoveParentEdge.parentNode;
            secondNodeToRemoveParent.removeChildEdge(secondNodeToRemoveParentEdge);
            secondNodeToRemove.removeParentEdge(secondNodeToRemoveParentEdge);
            secondNodeToRemoveParent.addChildEdge(secondEdgeToExtend);
        }

        return logHR;
    }


    /**
     * Find removable reassortment edges (single-segment, narrow topology)
     */
    private List<NetworkEdge> findRemovableNarrowReassortments() {
        List<NetworkEdge> removable = new ArrayList<>();

        for (NetworkEdge edge : networkEdges) {
            if (edge.isRootEdge()) continue;

            // Must carry exactly 1 segment
            if (edge.hasSegments.cardinality() != 1) continue;

            // Must be child of reassortment, parent of coalescence
            if (!edge.childNode.isReassortment()) continue;
            if (!edge.parentNode.isCoalescence()) continue;

            // Check narrow topology: the spouse edge should also connect to same coalescence
            NetworkEdge spouse = getSpouseEdge(edge);
            if (spouse.parentNode == edge.parentNode) {
                removable.add(edge);
            }
        }

        return removable;
    }

    /**
     * Count removable reassortments
     */
    private int countRemovableNarrowReassortments() {
        return findRemovableNarrowReassortments().size();
    }

    /**
     * Select a random segment from a BitSet
     */
    private int selectRandomSegment(BitSet segments) {
        int count = segments.cardinality();
        int selected = Randomizer.nextInt(count);
        int segment = -1;

        for (int i = 0; i <= selected; i++) {
            segment = segments.nextSetBit(segment + 1);
        }

        return segment;
    }

    /**
     * Check if node has at least one internal child (used for counting valid grandparents)
     * From NetworkExchange.isg()
     */
    private int isg(final Node n) {
        return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
    }

    /**
     * Find the network edge at a specific height for a given segment
     * From NetworkExchange.getEdgeWithHeight()
     */
    private NetworkEdge getEdgeWithHeight(NetworkEdge edge, double height, int segment) {
        if (edge.childNode.getHeight() <= height && edge.parentNode.getHeight() >= height) {
            return edge;
        } else {
            for (NetworkEdge parentEdge : edge.parentNode.getParentEdges()) {
                if (parentEdge.hasSegments.get(segment)) {
                    return getEdgeWithHeight(parentEdge, height, segment);
                }
            }
        }
        return null;
    }

    /**
     * Get the segment tree node number for a network edge
     * Traverses down the network edge to find the corresponding tree node
     */
    private int getSegmentTreeNodeNr(NetworkEdge edge, int segment) {
        NetworkNode node = edge.childNode;

        // Traverse down to find the segment tree node
        while (node != null) {
            if (node.isLeaf()) {
                return node.segmentIndices[segment];
            }

            if (node.isCoalescence()) {
                // Check if both children have the segment (coalescence for this segment)
                if (node.getChildEdges().get(0).hasSegments.get(segment) &&
                    node.getChildEdges().get(1).hasSegments.get(segment)) {
                    return node.segmentIndices[segment];
                }
                // Follow the branch that has the segment
                if (node.getChildEdges().get(0).hasSegments.get(segment)) {
                    node = node.getChildEdges().get(0).childNode;
                } else if (node.getChildEdges().get(1).hasSegments.get(segment)) {
                    node = node.getChildEdges().get(1).childNode;
                } else {
                    return -1;
                }
            } else if (node.isReassortment()) {
                // Follow down through reassortment
                node = node.getChildEdges().get(0).childNode;
            } else {
                return -1;
            }
        }

        return -1;
    }

    /**
     * Find network edge corresponding to a tree node number
     * Starts from a parent edge and searches downward for the tree node
     */
    private NetworkEdge findNetworkEdgeForTreeNode(NetworkEdge startEdge, int treeNodeNr, int segment) {
        // Recursively search down from startEdge to find edge with this tree node
        return findNetworkEdgeRecursive(startEdge.childNode, treeNodeNr, segment);
    }

    private NetworkEdge findNetworkEdgeRecursive(NetworkNode node, int treeNodeNr, int segment) {
        if (node.isLeaf()) {
            if (node.segmentIndices != null && node.segmentIndices[segment] == treeNodeNr) {
                // Return the edge above this leaf
                for (NetworkEdge e : networkEdges) {
                    if (e.childNode == node && e.hasSegments.get(segment)) {
                        return e;
                    }
                }
            }
            return null;
        }

        if (node.isCoalescence()) {
            // Check if this node corresponds to the tree node
            if (node.segmentIndices != null && node.segmentIndices[segment] == treeNodeNr) {
                // Check if both children have the segment (true coalescence for this segment)
                if (node.getChildEdges().get(0).hasSegments.get(segment) &&
                    node.getChildEdges().get(1).hasSegments.get(segment)) {
                    // Return edge above this node
                    for (NetworkEdge e : networkEdges) {
                        if (e.childNode == node && e.hasSegments.get(segment)) {
                            return e;
                        }
                    }
                }
            }

            // Search down both children
            NetworkEdge e0 = node.getChildEdges().get(0);
            NetworkEdge e1 = node.getChildEdges().get(1);

            if (e0.hasSegments.get(segment)) {
                NetworkEdge result = findNetworkEdgeRecursive(e0.childNode, treeNodeNr, segment);
                if (result != null) return result;
            }

            if (e1.hasSegments.get(segment)) {
                NetworkEdge result = findNetworkEdgeRecursive(e1.childNode, treeNodeNr, segment);
                if (result != null) return result;
            }

        } else if (node.isReassortment()) {
            // Follow down
            return findNetworkEdgeRecursive(node.getChildEdges().get(0).childNode, treeNodeNr, segment);
        }

        return null;
    }

}
