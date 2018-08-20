package coalre.operators;

import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class DivertSegmentOperator extends NetworkOperator {

    @Override
    public double networkProposal() {
        double logHR = 0.0;

        List<NetworkEdge> sourceEdges = network.getEdges().stream()
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.hasSegments.cardinality()>1)
                .collect(Collectors.toList());

        if (sourceEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        logHR -= Math.log(1.0/sourceEdges.size());

        NetworkEdge sourceEdge = sourceEdges.get(Randomizer.nextInt(sourceEdges.size()));
        NetworkEdge destEdge = getSpouseEdge(sourceEdge);

        BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
        logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments);


        NetworkNode[] mrcaNodes = getSegmentMRCAs(sourceEdge, destEdge, segsToDivert);

        network.startEditing(this);

        logHR += removeSegmentsFromAncestors(sourceEdge, segsToDivert, mrcaNodes);
        logHR -= addSegmentsToAncestors(destEdge, segsToDivert, mrcaNodes);

        if (!allEdgesAncestral())
            return Double.NEGATIVE_INFINITY;

        logHR += getLogConditionedSubsetProb(destEdge.hasSegments);

        int reverseSourceEdgeCount = (int)(network.getEdges().stream()
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.hasSegments.cardinality()>1)
                .count());

        logHR += Math.log(1.0/reverseSourceEdgeCount);

        return logHR;
    }


    /**
     * Remove segments from this edge and ancestors.
     *
     * @param edge edge at which to start removal
     * @param segsToRemove segments to remove from edge and ancestors
     * @return log probability of reverse operation
     */
    double removeSegmentsFromAncestors(NetworkEdge edge, BitSet segsToRemove, NetworkNode[] stopNodes) {
        double logP = 0.0;

        segsToRemove = (BitSet)segsToRemove.clone();
        segsToRemove.and(edge.hasSegments);

        if (segsToRemove.isEmpty())
            return logP;

        edge.hasSegments.andNot(segsToRemove);

        if (edge.isRootEdge())
            return logP;

        for (int segIdx=segsToRemove.nextSetBit(0); segIdx != -1;
             segIdx = segsToRemove.nextSetBit(segIdx)) {

            if (stopNodes[segIdx] != null && stopNodes[segIdx] == edge.parentNode)
                segsToRemove.clear(segIdx);
        }

        if (edge.parentNode.isReassortment()) {

            logP += Math.log(0.5)*segsToRemove.cardinality();

            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove, stopNodes);
            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(1), segsToRemove, stopNodes);

        } else {

            segsToRemove.andNot(getSisterEdge(edge).hasSegments);
            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove, stopNodes);

        }

        return logP;
    }

    /**
     * Add segments to this edge and ancestors.
     *
     * @param edge edge at which to start addition
     * @param segsToAdd segments to add to the edge and ancestors
     * @return log probability of operation
     */
    double addSegmentsToAncestors(NetworkEdge edge, BitSet segsToAdd, NetworkNode[] stopNodes) {
        double logP = 0.0;

        segsToAdd = (BitSet)segsToAdd.clone();
        segsToAdd.andNot(edge.hasSegments);

        if (segsToAdd.isEmpty())
            return logP;

        edge.hasSegments.or(segsToAdd);

        if (edge.isRootEdge())
            return logP;

        for (int segIdx = segsToAdd.nextSetBit(0); segIdx != -1; segIdx = segsToAdd.nextSetBit(segIdx)) {
            if (stopNodes[segIdx] != null && edge.parentNode == stopNodes[segIdx])
                segsToAdd.clear(segIdx);
        }

        if (edge.parentNode.isReassortment()) {

            BitSet segsToAddLeft = new BitSet();
            BitSet segsToAddRight = new BitSet();

            for (int segIdx=segsToAdd.nextSetBit(0); segIdx != -1;
                    segIdx=segsToAdd.nextSetBit(segIdx+1)) {
                if (Randomizer.nextBoolean())
                    segsToAddLeft.set(segIdx);
                else
                    segsToAddRight.set(segIdx);

                logP += Math.log(0.5);
            }

            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAddLeft, stopNodes);
            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(1), segsToAddRight, stopNodes);

        } else {

            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAdd, stopNodes);
        }

        return logP;
    }

    /**
     * Check that each edge is ancestral to at least one segment.
     *
     * @return true if all edges are ancestral.
     */
    public boolean allEdgesAncestral() {
        Set<NetworkNode> nodeList = networkInput.get().getNodes();
        for (NetworkNode node : nodeList) {
            for (NetworkEdge parentEdge : node.getParentEdges()) {
                if (parentEdge.hasSegments.isEmpty())
                    return false;
            }
        }

        return true;
    }

    public Set<NetworkNode> getAncestralNodes(NetworkEdge startEdge, int segIdx) {

        Set<NetworkNode> ancestralNodes = new HashSet<>();

        NetworkEdge thisEdge = startEdge;

        while (!thisEdge.isRootEdge()) {
            ancestralNodes.add(startEdge.parentNode);

            for (NetworkEdge parentEdge : thisEdge.parentNode.getParentEdges()) {
                if (parentEdge.hasSegments.get(segIdx)) {
                    thisEdge = parentEdge;
                    break;
                }
            }
        }

        return ancestralNodes;
    }

    public NetworkNode getSegmentMRCA(NetworkEdge edge1, NetworkEdge edge2, int segIdx) {

        Set<NetworkNode> ancestralNodes = getAncestralNodes(edge1, segIdx);

        NetworkEdge thisEdge = edge2;

        while (!thisEdge.isRootEdge()) {
            if (ancestralNodes.contains(thisEdge.parentNode))
                return thisEdge.parentNode;

            for (NetworkEdge parentEdge : thisEdge.parentNode.getParentEdges()) {
                if (parentEdge.hasSegments.get(segIdx)) {
                    thisEdge = parentEdge;
                    break;
                }
            }
        }

        return null;
    }

    public NetworkNode[] getSegmentMRCAs(NetworkEdge edge1, NetworkEdge edge2, BitSet segs) {
        NetworkNode[] mrcaNodes = new NetworkNode[network.getSegmentCount()];

        for (int segIdx=segs.nextSetBit(0); segIdx != -1; segIdx = segs.nextSetBit(segIdx)) {
            mrcaNodes[segIdx] = getSegmentMRCA(edge1, edge2, segIdx);
        }

        return mrcaNodes;
    }

}
