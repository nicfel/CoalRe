package coalre.operators;

import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

public class DivertSegmentOperator extends NetworkOperator {

    protected Network network;

    @Override
    public void initAndValidate() {
        network = networkInput.get();
    }

    @Override
    public double proposal() {
        double logHR = 0.0;

        List<NetworkNode> reassortmentNodes = network.getNodes().stream()
                .filter(NetworkNode::isReassortment)
                .filter(n -> n.getChildEdges().get(0).hasSegments.cardinality()>1)
                .collect(Collectors.toList());

        if (reassortmentNodes.isEmpty())
            return Double.NEGATIVE_INFINITY;

        NetworkNode reassortmentNode = reassortmentNodes.get(Randomizer.nextInt(reassortmentNodes.size()));
        logHR -= Math.log(1.0/reassortmentNodes.size());

        NetworkEdge sourceEdge = reassortmentNode.getParentEdges().get(Randomizer.nextInt(2));
        NetworkEdge siblingEdge = getSpouseEdge(sourceEdge);

        BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);

        return logHR;
    }


    /**
     * Remove segments from this edge and ancestors.
     *
     * @param edge edge at which to start removal
     * @param segsToRemove segments to remove from edge and ancestors
     * @return log probability of reverse operation
     */
    double removeSegmentsFromAncestors(NetworkEdge edge, BitSet segsToRemove) {
        double logP = 0.0;

        if (!edge.hasSegments.intersects(segsToRemove))
            return logP;

        segsToRemove = (BitSet)segsToRemove.clone();
        segsToRemove.and(edge.hasSegments);

        edge.hasSegments.andNot(segsToRemove);

        if (edge.isRootEdge())
            return logP;

        if (edge.parentNode.isCoalescence()) {
            segsToRemove.andNot(getSisterEdge(edge).hasSegments);
        }

        if (edge.parentNode.isReassortment())
            logP += -LOG2*segsToRemove.cardinality();

        for (NetworkEdge parentEdge : edge.parentNode.getParentEdges())
            logP += removeSegmentsFromAncestors(parentEdge, segsToRemove);

        return logP;
    }

    /**
     * Add segments to this edge and ancestors.
     *
     * @param edge edge at which to start addition
     * @param segsToAdd segments to add to the edge and ancestors
     * @return log probability of operation
     */
    double addSegmentsToAncestors(NetworkEdge edge, BitSet segsToAdd) {
        double logP = 0.0;

        segsToAdd = (BitSet)segsToAdd.clone();
        segsToAdd.andNot(edge.hasSegments);

        if (segsToAdd.isEmpty())
            return logP;

        edge.hasSegments.or(segsToAdd);

        if (edge.isRootEdge())
            return logP;

        if (edge.parentNode.isReassortment()) {

            BitSet segsToAddLeft = new BitSet();
            BitSet segsToAddRight = new BitSet();

            for (int segIdx=segsToAdd.nextSetBit(0); segIdx != -1;
                    segIdx=segsToAdd.nextSetBit(segIdx+1)) {
                if (Randomizer.nextBoolean())
                    segsToAddLeft.set(segIdx);
                else
                    segsToAddRight.set(segIdx);

                logP += -LOG2;
            }

            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAddLeft);
            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(1), segsToAddRight);

        } else {

            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAdd);
        }

        return logP;
    }

    /**
     * Check that each edge is ancestral to at least one segment.
     *
     * @return true if all edges are ancestral.
     */
    public boolean allEdgesAncestral() {
        for (NetworkNode node : networkInput.get().getNodes()) {
            for (NetworkEdge parentEdge : node.getParentEdges()) {
                if (parentEdge.hasSegments.isEmpty())
                    return false;
            }
        }

        return true;
    }

}
