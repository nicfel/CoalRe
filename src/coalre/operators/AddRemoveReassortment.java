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

            return addReassortment();

        } else {

            return removeReassortment();

        }
    }

    private double addReassortment() {

        return 0.0;
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


        return 0.0;
    }

    /**
     * Remove segments from this edge and ancestors.
     *
     * @param edge edge at which to start removal
     * @param segsToRemove segments to remove from edge and ancestors
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
     * Retrieve sister edge given parent.
     * @param childEdge child edge
     * @return sister of given child edge
     */
    private NetworkEdge getSisterEdge(NetworkEdge childEdge) {
        int idx = childEdge.parentNode.getChildEdges().indexOf(childEdge);
        int otherIdx = (idx + 1) % 2;

        return childEdge.parentNode.getChildEdges().get(otherIdx);
    }
}
