package coalre.operators;

import beast.core.Input;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Set;

public class AddRemoveReassortment extends NetworkOperator {

    public Input<Double> alphaInput = new Input<>("alpha",
            "Mean of exponential used for choosing root attachment times.",
            Input.Validate.REQUIRED);

    private Network network;
    private double alpha;

    @Override
    public void initAndValidate() {
        network = networkInput.get();
        alpha = alphaInput.get();
    }

    int count = 0;

    @Override
    public double proposal() {

        System.out.print("Add/remove counter " + (count++));


        if (Randomizer.nextBoolean()) {

            System.out.println(" (ADD)");

            return addReassortment();

        } else {

            System.out.println(" (REMOVE)");

            return removeReassortment();

        }
    }

    double addReassortment() {
        double logHR = 0.0;

        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        NetworkEdge sourceEdge;
        do {
            sourceEdge = networkEdges.get(Randomizer.nextInt(networkEdges.size()));
        } while (sourceEdge.isRootEdge());

        if (sourceEdge.hasSegments.cardinality()<2) {
            System.out.println("DIRECT REJECT (TOO FEW SEGS ON SELECTED EDGE)");
            return Double.NEGATIVE_INFINITY;
        }

        double sourceTime = Randomizer.nextDouble()*sourceEdge.getLength() + sourceEdge.childNode.getHeight();

        logHR -= Math.log(1.0/(networkEdges.size()-1)/sourceEdge.getLength());

        NetworkEdge destEdge = networkEdges.get(Randomizer.nextInt(networkEdges.size()));

        if (!destEdge.isRootEdge() && sourceTime>destEdge.parentNode.getHeight()) {
            System.out.println("DIRECT REJECT (INCOMPATIBLE SOURCE/DEST PAIR)");
            return Double.NEGATIVE_INFINITY;
        }

        logHR -= Math.log(1.0/networkEdges.size());

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        double destTime;
        if (destEdge.isRootEdge()) {

            destTime = minDestTime + Randomizer.nextExponential(1.0/alpha);
            logHR -= -(1.0/alpha)*(destTime-minDestTime) + Math.log(1.0/alpha);

        } else {

            destTime = Randomizer.nextDouble()*(destEdge.parentNode.getHeight()-minDestTime) + minDestTime;
            logHR -= Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));

        }

        network.startEditing(this);

        // Create new reassortment edge

        NetworkNode sourceNode = new NetworkNode();
        sourceNode.setHeight(sourceTime);

        NetworkNode oldSourceEdgeParent = sourceEdge.parentNode;
        oldSourceEdgeParent.removeChildEdge(sourceEdge);
        sourceNode.addChildEdge(sourceEdge);

        NetworkEdge newEdge1 = new NetworkEdge();
        sourceNode.addParentEdge(newEdge1);
        oldSourceEdgeParent.addChildEdge(newEdge1);

        newEdge1.hasSegments = (BitSet) sourceEdge.hasSegments.clone();

        if (destEdge == sourceEdge)
            destEdge = newEdge1;

        NetworkNode destNode = new NetworkNode();
        destNode.setHeight(destTime);

        NetworkNode oldDestEdgeParent = destEdge.parentNode;
        if (oldDestEdgeParent != null) {
            oldDestEdgeParent.removeChildEdge(destEdge);
        }

        destNode.addChildEdge(destEdge);

        NetworkEdge newEdge2 = new NetworkEdge();
        destNode.addParentEdge(newEdge2);

        if (oldDestEdgeParent == null) {
            network.setRootEdge(newEdge2);
        } else {
            oldDestEdgeParent.addChildEdge(newEdge2);
        }

        newEdge2.hasSegments = (BitSet) destEdge.hasSegments.clone();

        NetworkEdge reassortmentEdge = new NetworkEdge();
        sourceNode.addParentEdge(reassortmentEdge);
        destNode.addChildEdge(reassortmentEdge);
        reassortmentEdge.hasSegments = new BitSet();

        // Choose segments to divert to new edge
        BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
        logHR += removeSegmentsFromAncestors(newEdge1, segsToDivert);
        logHR -= addSegmentsToAncestors(reassortmentEdge, segsToDivert);

        if (!allEdgesAncestral()) {
            System.out.println("DIRECT REJECT (EDGE WITH NO SEGMENT)");
            return Double.NEGATIVE_INFINITY;
        }

        // HR contribution for reverse move
        int nRemovableEdges = (int) network.getEdges().stream()
                .filter(e -> !e.isRootEdge()
                        && e.childNode.isReassortment()
                        && e.parentNode.isCoalescence()).count();
        logHR += Math.log(1.0/nRemovableEdges);

        return logHR;
    }

    double removeReassortment() {
        double logHR = 0.0;

        Set<NetworkEdge> networkEdges = network.getEdges();
        List<NetworkEdge> removableEdges = new ArrayList<>();

        for (NetworkEdge edge : networkEdges)
            if (!edge.isRootEdge() && edge.childNode.isReassortment() && edge.parentNode.isCoalescence())
                removableEdges.add(edge);

        if (removableEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        network.startEditing(this);

        NetworkEdge edgeToRemove = removableEdges.get(Randomizer.nextInt(removableEdges.size()));
        NetworkNode nodeToRemove = edgeToRemove.childNode;
        NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
        NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;

        logHR -= Math.log(1.0/(removableEdges.size()));

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

        if (!allEdgesAncestral()) {
            System.out.println("DIRECT REJECT (EDGE WITH NO SEGMENT)");
            return Double.NEGATIVE_INFINITY;
        }

        // HR contribution of edge selection for reverse move
        logHR += Math.log(1.0/(network.getEdges().size()-1)/edgeToExtend.getLength());

        // TODO: HR contributions for time selections in reverse move.

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
        for (NetworkNode node : network.getNodes()) {
            for (NetworkEdge parentEdge : node.getParentEdges()) {
                if (parentEdge.hasSegments.isEmpty())
                    return false;
            }
        }

        return true;
    }

}
