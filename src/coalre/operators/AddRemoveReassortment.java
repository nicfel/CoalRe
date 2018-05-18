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

    private Network network;

    @Override
    public void initAndValidate() {
        network = networkInput.get();
    }

    @Override
    public double proposal() {


        if (Randomizer.nextBoolean()) {

            return addReassortment();

        } else {

            return removeReassortment();

        }
    }

    double addReassortment() {

        return 0.0;
    }

    double removeReassortment() {
        double logHR = 0.0;

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
        NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
        NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;

        logHR -= Math.log(1.0/(2*reassortmentNodes.size()));

        // Divert segments away from chosen edge
        BitSet segsToDivert = (BitSet)edgeToRemove.hasSegments.clone();
        logHR += removeSegmentsFromAncestors(edgeToRemove, segsToDivert);
        logHR -= addSegmentsToAncestors(edgeToRemoveSpouse, segsToDivert);

        // Remove edge and associated nodes
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

        if (secondNodeToRemove.getParentEdges().get(0).isRootEdge()) {
            network.setRootEdge(secondEdgeToExtend);

        } else {
            NetworkEdge secondNodeToRemoveParentEdge = secondNodeToRemove.getParentEdges().get(0);
            NetworkNode secondNodeToRemoveParent = secondNodeToRemoveParentEdge.parentNode;
            secondNodeToRemoveParent.removeChildEdge(secondNodeToRemoveParentEdge);
            secondNodeToRemove.removeParentEdge(secondNodeToRemoveParentEdge);

            secondNodeToRemoveParent.addChildEdge(secondEdgeToExtend);
        }

        // HR contribution of edge selection for reverse move
        logHR += Math.log(1.0/(network.getNodes().size()-1));

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

}
