package coalre.operators;

import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Set;

public class AddRemoveReassortment extends NetworkOperator {

    Network network;

    @Override
    public void initAndValidate() {
        network = networkInput.get();
    }

    @Override
    public double proposal() {


        if (Randomizer.nextBoolean()) {
            return addReassortment(List<NetworkNode>);

        } else {
            return removeReassortment();
        }

        return 0;
    }

    private double addReassortment() {

    }

    private double removeReassortment() {
        Set<NetworkNode> networkNodes = network.getNodes();
        List<NetworkNode> reassortmentNodes = new ArrayList<>();
        for (NetworkNode node : networkNodes)
            if (node.isReassortment())
                reassortmentNodes.add(node);

        if (reassortmentNodes.isEmpty())
            return Double.NEGATIVE_INFINITY;

        NetworkNode nodeToRemove = reassortmentNodes.get(Randomizer.nextInt(reassortmentNodes.size()));
        int removalEdgeIdx = Randomizer.nextInt(2);
        NetworkEdge edgeToRemove = nodeToRemove.getParentEdges().get(removalEdgeIdx);

        // Remove edge parent node
        edgeToRemove.parentNode

    }

    /**
     * Remove segments from this edge and ancestors.
     *
     * @param edge edge at which to start removal
     * @param segsToRemove segments to remove from edge and ancestors
     * @param seenNodes set of nodes already seen during traversal
     */
    private void removeSegmentsFromAncestors(NetworkEdge edge, BitSet segsToRemove, Set<NetworkNode> seenNodes) {
        if (!edge.hasSegments.intersects(segsToRemove))
            return;

        edge.hasSegments.andNot(segsToRemove);

        if (edge.parentNode.isCoalescence()) {
            segsToRemove = (BitSet)segsToRemove.clone();
            segsToRemove.andNot(getSisterEdge(edge.parentNode, edge).hasSegments);
        }

         if (edge.isRootEdge() || seenNodes.contains(edge.parentNode))
            return;

         seenNodes.add(edge.parentNode);

        for (NetworkEdge parentEdge : edge.parentNode.getParentEdges())
            removeSegmentsFromAncestors(parentEdge, segsToRemove, seenNodes);
    }

    private void addSegmentsToAncestors(NetworkEdge edge, BitSet segsToAdd, Set<NetworkNode> seenNodes) {
        int origSegCount = edge.hasSegments.cardinality();
        edge.hasSegments.or(segsToAdd);

        // Stop here if nothing has changed
        if (edge.hasSegments.cardinality() == origSegCount)
            return;

        if (edge.isRootEdge() || seenNodes.contains(edge.parentNode))
            return;


        seenNodes.add(edge.parentNode);

        BitSet segsToAddLeft = new BitSet();
        BitSet segsToAddRight = new BitSet();
    }

    /**
     * Retrieve sister edge given parent.
     * @param parentNode parent node
     * @param childEdge child edge
     * @return sister of given child edge
     */
    private NetworkEdge getSisterEdge(NetworkNode parentNode, NetworkEdge childEdge) {
        int idx = (parentNode.getChildEdges().indexOf(childEdge);
        int otherIdx = (idx + 1) % 2;

        return parentNode.getChildEdges().get(otherIdx);
    }
}
