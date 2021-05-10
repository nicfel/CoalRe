package coalre.operators;

import beast.core.Description;
import beast.core.Input;
import beast.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.*;
import java.util.stream.Collectors;
@Description("Does add and remove reassortment events by sampling the reattachment heights from and exponential distribution." +
			"The operation is tunable. In the future, it could use a coalescent style way to reattach reassortment nodes")
public class LocalAddRemoveReassortment extends DivertSegmentOperator {

    public Input<Double> alphaInput = new Input<>("alpha",
            "Mean of exponential used for choosing root attachment times.",
            Input.Validate.REQUIRED);

    private double alpha;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        alpha = alphaInput.get();
    }

    @Override
    public double networkProposal() {

        double logHR;
        if (Randomizer.nextBoolean())
            logHR = addReassortment();
        else
            logHR = removeReassortment();

        return logHR;
    }

    double addReassortment() {
        double logHR = 0.0;

        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        List<NetworkEdge> possibleSourceEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=2)
                .collect(Collectors.toList());

        NetworkEdge sourceEdge = possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size()));
        double sourceTime = Randomizer.nextDouble()*sourceEdge.getLength() + sourceEdge.childNode.getHeight();

        logHR -= Math.log(1.0/(double) possibleSourceEdges.size())
                + Math.log(1.0/sourceEdge.getLength());
        
        
        // randomly sample the destination time from and exponential
        double destTime = sourceTime + Randomizer.nextExponential(1.0/alpha);
        
        // add the height HR contribution
        logHR -= -(1.0/alpha)*(destTime-sourceTime) + Math.log(1.0/alpha);

        
//        List<NetworkEdge> possibleDestinationEdges = new ArrayList<>();
//        // check if the destination time is greater that the root edge time
//        if (destTime>network.getRootEdge().childNode.getHeight()){
//        	possibleDestinationEdges.add(network.getRootEdge());
//        }else{        
//	        // check which are possible destination edges there are with that height
//        	possibleDestinationEdges = networkEdges.stream()
//	                .filter(e -> !e.isRootEdge())
//	                .filter(e -> e.childNode.getHeight()<destTime)
//	                .filter(e -> e.parentNode.getHeight()>destTime)
//	                .collect(Collectors.toList());
//        }

//        NetworkEdge destEdge = possibleDestinationEdges.get(Randomizer.nextInt(possibleDestinationEdges.size()));
        
        NetworkEdge destEdge = networkEdges.get(Randomizer.nextInt(networkEdges.size()));
        logHR -= Math.log(1.0/networkEdges.size());
        
        if (destEdge.childNode.getHeight()>destTime)
        	return Double.NEGATIVE_INFINITY;
        
        if (!destEdge.isRootEdge() && destEdge.parentNode.getHeight()<destTime)
        	return Double.NEGATIVE_INFINITY;

        
        
//        logHR -= Math.log(1.0/possibleDestinationEdges.size());
        

        logHR += addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // HR contribution for reverse move
        int nRemovableEdges = (int) network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=1)
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.parentNode.isCoalescence())
                .count();
        
        logHR += Math.log(1.0/nRemovableEdges);

        return logHR;
    }

    double addReassortmentEdge(NetworkEdge sourceEdge, double sourceTime,
                               NetworkEdge destEdge, double destTime) {

        double logHR = 0.0;

        network.startEditing(this);

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
        logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments);
        logHR -= addSegmentsToAncestors(reassortmentEdge, segsToDivert);
        logHR += removeSegmentsFromAncestors(newEdge1, segsToDivert);

        return logHR;
    }

    double removeReassortment() {
        double logHR = 0.0;

        List<NetworkEdge> removableEdges = network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=1)
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.parentNode.isCoalescence())
                .collect(Collectors.toList());

        if (removableEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        // pick a random edge to remove
        NetworkEdge edgeToRemove = removableEdges.get(Randomizer.nextInt(removableEdges.size()));
        logHR -= Math.log(1.0/(removableEdges.size()));

        double sourceTime = edgeToRemove.childNode.getHeight();
        
        NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
        
        NetworkEdge destEdge = getSisterEdge(edgeToRemove);
        if (destEdge.childNode == edgeToRemove.childNode)
            destEdge = sourceEdge;
        double destTime = edgeToRemove.parentNode.getHeight();
               


        // Remove reassortment edge
        logHR += removeReassortmentEdge(edgeToRemove);

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // HR contribution for reverse
        int nPossibleSourceEdges = (int) network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=2)
                .count();

        // reverse HR contribution for the edge length of the source edge
        logHR += Math.log(1.0/(double) nPossibleSourceEdges)
                + Math.log(1.0/sourceEdge.getLength());
        
        // reverse HR contribution for height reattachment
        logHR += -(1.0/alpha)*(destTime-sourceTime) + Math.log(1.0/alpha);

        logHR += Math.log(1.0/network.getEdges().size());

        
//        // check how many possible destinations there are
//        int nPossibleDestEdges = 1;
//        if (destTime<network.getRootEdge().childNode.getHeight())
//        	nPossibleDestEdges = (int) network.getEdges().stream()
//                .filter(e -> !e.isRootEdge())
//                .filter(e -> e.childNode.getHeight()<destTime)
//                .filter(e -> e.parentNode.getHeight()>destTime)
//                .count();
//
//        // reverse HR contribution for the number of possible destination
//        logHR += Math.log(1.0/nPossibleDestEdges);

        return logHR;
    }


    double removeReassortmentEdge(NetworkEdge edgeToRemove) {
        double logHR = 0.0;

        network.startEditing(this);

        NetworkNode nodeToRemove = edgeToRemove.childNode;
        NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
        NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;

        // Divert segments away from chosen edge
        BitSet segsToDivert = (BitSet) edgeToRemove.hasSegments.clone();
        logHR -= addSegmentsToAncestors(edgeToRemoveSpouse, segsToDivert);
        logHR += removeSegmentsFromAncestors(edgeToRemove, segsToDivert);
        logHR += getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments);

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

        if (!networkTerminatesAtMRCA())
            return Double.NEGATIVE_INFINITY;

        return logHR;
    }

    /**
     * Simple (but probably too expensive) check for a kind of invalid network
     * which can result from an edge deletion operation: one in which the
     * network posesses nontrivial structure above the MRCA. (I.e. the MRCA
     * is not the root.)
     *
     * @return true if the network terminates at the true MRCA. (I.e. is valid.)
     */
    protected boolean networkTerminatesAtMRCA() {
        List<NetworkNode> sortedNodes = new ArrayList<>(network.getNodes());
        sortedNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));
        List<NetworkNode> sampleNodes = sortedNodes.stream().filter(NetworkNode::isLeaf).collect(Collectors.toList());
        double maxSampleHeight = sampleNodes.get(sampleNodes.size()-1).getHeight();

        int lineages = 0;
        for (NetworkNode node : sortedNodes) {
            switch(node.getChildEdges().size()) {
                case 2:
                    // Coalescence

                    lineages -= 1;
                    break;

                case 1:
                    // Reassortment

                    if (lineages < 2 && node.getHeight() > maxSampleHeight)
                        return false;

                    lineages += 1;
                    break;

                case 0:
                    // Sample

                    lineages += 1;
                    break;
            }
        }

        return true;
    }
}
