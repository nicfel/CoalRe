package coalre.operators;

import beast.core.Input;
import beast.util.Package;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

public class AddRemoveWithEmptySegment extends DivertSegmentOperator {

    public Input<Double> alphaInput = new Input<>("alpha",
            "Mean of exponential used for choosing root attachment times.",
            Input.Validate.REQUIRED);
    
    public Input<Double> lambdaInput = new Input<>("lambda",
            "lambda of th epoisson distribution for how many empty edges to add.",
            50.0);


    private double alpha;
    private double lambda;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        alpha = alphaInput.get();
        lambda = lambdaInput.get();
    }

    @Override
    public double networkProposal() {
//    	System.out.println("new proposal");

        double logHR = 0.0;
        
        // Adds empty network edges (for reversibly removing all empty network edges)
//        System.out.println(network.getExtendedNewick());

        addEmptyNetworkSegments();
        
//        System.out.println(network.getExtendedNewick());
        
        
        if (Randomizer.nextBoolean())
            logHR += addReassortment();
        else
            logHR += removeReassortment();
        
//        System.out.println(network.getExtendedNewick());

        // removes all the empty network segments again
//        System.out.println(network.getExtendedNewick());
        RemoveAllEmptyNetworkSegments();
//        System.out.println(network.getExtendedNewick());
//        System.exit(0);

        return logHR;
    }
    
    private double addEmptyNetworkSegments(){
    	double logHR = 0.0;
    	// randomly sample the number of edges to add
    	int nrEmptyEdges = (int) Randomizer.nextPoisson(lambdaInput.get());
    	// TODO prob contribution of that
    	for (int i = 0; i < nrEmptyEdges; i ++){
    		logHR += addEmptyReassortment();
    	}
    	return logHR;
    }
    
    
    double addEmptyReassortment() {
        double logHR = 0.0;

        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        List<NetworkEdge> possibleSourceEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=1)
                .collect(Collectors.toList());

        NetworkEdge sourceEdge = possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size()));
        double sourceTime = Randomizer.nextDouble()*sourceEdge.getLength() + sourceEdge.childNode.getHeight();

        logHR -= Math.log(1.0/(double)possibleSourceEdges.size())
                + Math.log(1.0/sourceEdge.getLength());

        NetworkEdge destEdge = networkEdges.get(Randomizer.nextInt(networkEdges.size()));
        // TODO make that proper
        int count=0;
        while (destEdge.hasSegments.cardinality()>=1){
        	destEdge = networkEdges.get(Randomizer.nextInt(networkEdges.size()));
        	count++;
        	if (count>5)
        		return Double.NEGATIVE_INFINITY;
        }
        
        logHR -= Math.log(1.0/networkEdges.size());

        if (!destEdge.isRootEdge() && destEdge.parentNode.getHeight() < sourceTime)
            return Double.NEGATIVE_INFINITY;

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        double destTime;
        if (destEdge.isRootEdge()) {

            destTime = minDestTime + Randomizer.nextExponential(1.0/alpha);
            logHR -= -(1.0/alpha)*(destTime-minDestTime) + Math.log(1.0/alpha);

        } else {

            destTime = Randomizer.nextDouble()*(destEdge.parentNode.getHeight()-minDestTime) + minDestTime;
            logHR -= Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));

        }

        // Create new reassortment edge

        logHR += addEmptyReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // HR contribution for reverse move
        int nRemovableEdges = (int) network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.parentNode.isCoalescence())
                .count();
        logHR += Math.log(1.0/nRemovableEdges);

        return logHR;
    }
    
    // Only adds the reassortment edge, but does not diverge segments
    double addEmptyReassortmentEdge(NetworkEdge sourceEdge, double sourceTime,
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
//		BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
//		logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments);
//		logHR -= addSegmentsToAncestors(reassortmentEdge, segsToDivert);
//		logHR += removeSegmentsFromAncestors(newEdge1, segsToDivert);
		
		//if (!allEdgesAncestral()){
		////return Double.NEGATIVE_INFINITY;
		//logHR += removeAllEmptyEdgesAncestral();
		//}
		
		return logHR;
	}

    
    public double RemoveAllEmptyNetworkSegments() {
    	double logHr = 0.0;
    	
        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        List<NetworkEdge> edgesToRemove = networkEdges.stream()
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.hasSegments.cardinality()==0)
                .collect(Collectors.toList());
        
        // TODO, see if ordering is costly, if yes, see if it's needed
        while (edgesToRemove.size()>0){
        	// get the node with the maximal height
        	int max_node = 0;
        	double max_height = 0.0;
        	for (int i = 0; i < edgesToRemove.size(); i++){
        		if (edgesToRemove.get(i).childNode.getHeight()>max_height){
        			max_height = edgesToRemove.get(i).childNode.getHeight();
        			max_node = i;
        		}
        	}
        		
        	removeEmptyReassortmentEdge(edgesToRemove.get(max_node));
            networkEdges = new ArrayList<>(network.getEdges());

            edgesToRemove = networkEdges.stream()
                    .filter(e -> e.childNode.isReassortment())
                    .filter(e -> e.hasSegments.cardinality()==0)
                    .collect(Collectors.toList());
        }
        
        if (!allEdgesAncestral()){
        	System.out.println("still has empty segments");
        	System.exit(0);
        	return Double.NEGATIVE_INFINITY;
        }

    	
//    	System.out.println("........................");
    	// goes through the network and removes all empty edges.
//        Set<NetworkNode> nodeList = networkInput.get().getNodes();
//        for (NetworkNode node : nodeList) {
//            for (NetworkEdge parentEdge : node.getParentEdges()) {
//                if (parentEdge.hasSegments.isEmpty()){
////                	System.out.println("removed something");
//                	System.out.println(network.getExtendedNewick());
//                	logHr -= removeEmptyReassortmentEdge(parentEdge);     
////                	System.out.println(node.getParentEdges());
//                	if (node.getParentEdges().size()==0)
//                		break;
//                }
//            }
//        }
        return logHr;
    }
    
    double removeEmptyReassortmentEdge(NetworkEdge edgeToRemove) {
        double logHR = 0.0;

//        NetworkEdge edgeToRemove = edgeToRemove;
//        logHR -= Math.log(1.0/(removableEdges.size()));

        double sourceTime = edgeToRemove.childNode.getHeight();
        NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
        NetworkEdge destEdge = getSisterEdge(edgeToRemove);
        if (destEdge.childNode == edgeToRemove.childNode)
            destEdge = sourceEdge;
        double destTime = edgeToRemove.parentNode.getHeight();

        // Remove reassortment edge
        logHR += actuallyRemoveEmptyReassortmentEdge(edgeToRemove);

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // HR contribution for reverse move

        Set<NetworkEdge> finalNetworkEdges = network.getEdges();

        int nPossibleSourceEdges = (int)finalNetworkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=2)
                .count();

        logHR += Math.log(1.0/(double)nPossibleSourceEdges)
                + Math.log(1.0/sourceEdge.getLength());

        logHR += Math.log(1.0/finalNetworkEdges.size());

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        if (destEdge.isRootEdge()) {
            logHR += -(1.0/alpha)*(destTime-minDestTime) + Math.log(1.0/alpha);
        } else {
            logHR += Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));
        }

        return logHR;
    }


    double actuallyRemoveEmptyReassortmentEdge(NetworkEdge edgeToRemove) {

        double logHR = 0.0;

        network.startEditing(this);

        NetworkNode nodeToRemove = edgeToRemove.childNode;
        NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
        NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;

//        // Divert segments away from chosen edge
//        BitSet segsToDivert = (BitSet) edgeToRemove.hasSegments.clone();
//        logHR -= addSegmentsToAncestors(edgeToRemoveSpouse, segsToDivert);
//        logHR += removeSegmentsFromAncestors(edgeToRemove, segsToDivert);
//        logHR += getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments);

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

//        if (!allEdgesAncestral()){
//            return Double.NEGATIVE_INFINITY;
//        }
        if (!networkTerminatesAtMRCA()){
            return Double.NEGATIVE_INFINITY;
        }



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

        logHR -= Math.log(1.0/(double)possibleSourceEdges.size())
                + Math.log(1.0/sourceEdge.getLength());

        NetworkEdge destEdge = networkEdges.get(Randomizer.nextInt(networkEdges.size()));
        logHR -= Math.log(1.0/networkEdges.size());

        if (!destEdge.isRootEdge() && destEdge.parentNode.getHeight() < sourceTime)
            return Double.NEGATIVE_INFINITY;

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        double destTime;
        if (destEdge.isRootEdge()) {

            destTime = minDestTime + Randomizer.nextExponential(1.0/alpha);
            logHR -= -(1.0/alpha)*(destTime-minDestTime) + Math.log(1.0/alpha);

        } else {

            destTime = Randomizer.nextDouble()*(destEdge.parentNode.getHeight()-minDestTime) + minDestTime;
            logHR -= Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));

        }

        // Create new reassortment edge

        logHR += addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // HR contribution for reverse move
        int nRemovableEdges = (int) network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
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

//        if (!allEdgesAncestral())
//            return Double.NEGATIVE_INFINITY;

        return logHR;
    }

    double removeReassortment() {
        double logHR = 0.0;

        List<NetworkEdge> removableEdges = network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.parentNode.isCoalescence())
                .collect(Collectors.toList());

        if (removableEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

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

        // HR contribution for reverse move

        Set<NetworkEdge> finalNetworkEdges = network.getEdges();

        int nPossibleSourceEdges = (int)finalNetworkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=2)
                .count();

        logHR += Math.log(1.0/(double)nPossibleSourceEdges)
                + Math.log(1.0/sourceEdge.getLength());

        logHR += Math.log(1.0/finalNetworkEdges.size());

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        if (destEdge.isRootEdge()) {
            logHR += -(1.0/alpha)*(destTime-minDestTime) + Math.log(1.0/alpha);
        } else {
            logHR += Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));
        }

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

        if (!networkTerminatesAtMRCA()){
            return Double.NEGATIVE_INFINITY;
        }



        return logHR;
    }
    
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
