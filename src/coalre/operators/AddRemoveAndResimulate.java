package coalre.operators;

import beast.base.core.Input;
import beast.base.util.Randomizer;
import coalre.distribution.CoalescentWithReassortment;
import coalre.distribution.NetworkEvent;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import java.util.*;
import java.util.stream.Collectors;

public class AddRemoveAndResimulate extends DivertSegmentAndResimulate {
    
	public Input<Boolean> randomlySampleAttachmentEdgeInput = new Input<>("randomlySampleAttachmentEdge",
			"Randomly sample edge to attach to", true);
		
    CoalescentWithReassortment coalescentDistr;
    
    double addProb = 0.5;
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
        coalescentDistr = coalescentDistrInput.get();
    }

    @Override
    public double networkProposal() {
        double logHR;
        
        for (int i = 0; i < segmentTrees.size(); i++) {
    		segmentsChanged.set(i, false);
        }


        if (Randomizer.nextDouble()< addProb) {
    		logHR = addRecombination();        	
        	logHR += Math.log((1.0 - addProb)/addProb);
        }else {
			logHR = removeRecombination();
            logHR += Math.log(addProb/(1.0 - addProb));
        }
        
        return logHR;
    }

    double addRecombination() {
        double logHR = 0.0;

        List<Integer> possibleSourceEdges = new ArrayList<>();
		for (int i = 0; i < networkEdges.size(); i++) {
			NetworkEdge edge = networkEdges.get(i);
			if (!edge.isRootEdge() && edge.hasSegments.cardinality() >= 2) {
				possibleSourceEdges.add(i);
			}
		}
        
        NetworkEdge sourceEdge = null;
        double sourceTime = -1.0;
        
		if (randomlySampleAttachmentEdgeInput.get()) {
	        
	        sourceEdge = networkEdges.get(possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size())));
	        sourceTime = Randomizer.nextDouble() * sourceEdge.getLength() + sourceEdge.childNode.getHeight();
	
	        logHR -= Math.log(1.0/(double)possibleSourceEdges.size())
	                + Math.log(1.0/sourceEdge.getLength());
		}else {
	        // compute the sum over all possible source edges
			double sumEdgeLengths = 0.0;
			for (Integer i : possibleSourceEdges) {
				sumEdgeLengths += networkEdges.get(i).getLength();
			}
					
			// Pick a random number between 0 and sum of edge lengths
			double randomEdgeLength = Randomizer.nextDouble() * sumEdgeLengths;
			logHR -= Math.log(1.0/sumEdgeLengths);
			// pick the source edge as the first edge whose length is greater than the random number
			double passedLength = 0;
			for (Integer i : possibleSourceEdges) {
				passedLength += networkEdges.get(i).getLength();
				if (passedLength > randomEdgeLength) {
					sourceEdge = networkEdges.get(i);
					sourceTime = passedLength-randomEdgeLength + sourceEdge.childNode.getHeight();
					break;
				}
			}

		}
                
    	// Calculate tree intervals
    	List<NetworkEvent> networkEventList = coalescentDistr.intervals.getNetworkEventList();

    	double currTime = sourceTime;
    	NetworkEvent prevEvent = null;
    	double attachmentTime = 0.0;
    	
    	for (NetworkEvent event : networkEventList) {
    		if (event.time>currTime) {
    			double rate = 0.5*prevEvent.lineages;
                double currentTransformedTime = coalescentDistr.populationFunction.getIntensity(currTime);
                double transformedTimeToNextCoal = Randomizer.nextExponential(rate);
                double timeToNextCoal = coalescentDistr.populationFunction.getInverseIntensity(
                        transformedTimeToNextCoal + currentTransformedTime);
                

                attachmentTime = timeToNextCoal;
                if (timeToNextCoal < event.time) {
                	logHR -= -rate * coalescentDistr.populationFunction.getIntegral(currTime, attachmentTime) +
                			Math.log(rate/coalescentDistr.populationFunction.getPopSize(attachmentTime));
                	break;
                }
                logHR -= -rate * coalescentDistr.populationFunction.getIntegral(currTime, event.time);
    			currTime = event.time;
    		}
    		prevEvent = event;
        }    	
    	
    	if (attachmentTime>network.getRootEdge().childNode.getHeight()) {    		
            double currentTransformedTime = coalescentDistr.populationFunction.getIntensity(network.getRootEdge().childNode.getHeight());
    		double transformedTimeToNextCoal = Randomizer.nextExponential(0.5);
            attachmentTime = coalescentDistr.populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime);
            
            logHR -= -0.5 * coalescentDistr.populationFunction.getIntegral(network.getRootEdge().childNode.getHeight(), attachmentTime);
        	logHR -= Math.log(0.5/coalescentDistr.populationFunction.getPopSize(attachmentTime));

    	}    	

    	double destTime = attachmentTime;
        // keep only those that coexist at the time of attachment
        List<NetworkEdge> destEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.parentNode.getHeight()>destTime)
                .filter(e -> e.childNode.getHeight()<=destTime)
               .collect(Collectors.toList());
        
        if (destEdges.size()==0)
        	destEdges.add(network.getRootEdge());
        
        
        NetworkEdge destEdge = destEdges.get(Randomizer.nextInt(destEdges.size()));
        logHR -= Math.log(1.0/destEdges.size());

        if (!destEdge.isRootEdge() && destEdge.parentNode.getHeight() < sourceTime)
            return Double.NEGATIVE_INFINITY;
      
        // Create new reassortment edge
        logHR += addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);
       
        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;  
        
        int nRemovableEdges = (int) networkEdges.stream().filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() >= 1).filter(e -> e.childNode.isReassortment())
				.filter(e -> e.parentNode.isCoalescence()).count();
        
        logHR += Math.log(1.0/nRemovableEdges);
        return logHR;
    }

    double addReassortmentEdge(NetworkEdge sourceEdge, double sourceTime,
                               NetworkEdge destEdge, double destTime) {

        double logHR = 0.0;

        boolean sourceisroot = sourceEdge.isRootEdge();
        NetworkNode sourceNode = new NetworkNode();
        sourceNode.setHeight(sourceTime);

        NetworkNode oldSourceEdgeParent = sourceEdge.parentNode;
        if (!sourceisroot)
        	oldSourceEdgeParent.removeChildEdge(sourceEdge);
        sourceNode.addChildEdge(sourceEdge);

        NetworkEdge newEdge1 = new NetworkEdge();
        sourceNode.addParentEdge(newEdge1);
        if (!sourceisroot)
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
        BitSet segsToDivert;
        if (divertOneSegmentInput.get()) {
        	double prob = 1./sourceEdge.hasSegments.cardinality();
            segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments, prob);
            logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments, segsToDivert, prob);
        }else {
	        segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
            logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments, segsToDivert, 0.5);
        }
        
        // Use resimulation approach instead of simple divertSegments
        logHR += divertSegmentsWithResimulation(reassortmentEdge, newEdge1, segsToDivert);
        
        networkEdges.add(reassortmentEdge);
        networkEdges.add(newEdge1);
        networkEdges.add(newEdge2);

        return logHR;
    }

    double removeRecombination() {
        double logHR = 0.0;
        
        List<Integer> removableEdges = new ArrayList<>();
		for (int i = 0; i < networkEdges.size(); i++) {
			NetworkEdge edge = networkEdges.get(i);
			if (!edge.isRootEdge()&& edge.childNode.isReassortment()
					&& edge.parentNode.isCoalescence() && edge.hasSegments.cardinality() >= 1 ) {
				removableEdges.add(i);
			}
		}
        
        if (removableEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        NetworkEdge edgeToRemove = networkEdges.get(removableEdges.get(Randomizer.nextInt(removableEdges.size())));
        logHR -= Math.log(1.0/(removableEdges.size()));
        

        double sourceTime = edgeToRemove.childNode.getHeight();
        NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
        NetworkEdge destEdge = getSisterEdge(edgeToRemove);
        if (destEdge.childNode == edgeToRemove.childNode)
            destEdge = sourceEdge;
        double destTime = edgeToRemove.parentNode.getHeight();        
    	// Calculate tree intervals
    	List<NetworkEvent> networkEventList = coalescentDistr.intervals.getNetworkEventList();

    	double currTime = sourceTime;
    	NetworkEvent prevEvent = null;
    	
    	for (NetworkEvent event : networkEventList) {
    		if (event.time>currTime) {                
            	double rate = 0.5 * (prevEvent.lineages-1);
                if (destTime<=event.time) {
                	logHR += -rate* coalescentDistr.populationFunction.getIntegral(currTime, destTime);                	
                	logHR += Math.log(rate/coalescentDistr.populationFunction.getPopSize(destTime));
                	break;
                }
                logHR += -rate * coalescentDistr.populationFunction.getIntegral(currTime, event.time);
    			currTime = event.time;
    		}
    		prevEvent = event;
        }    	
    	
        // Remove reassortment edge
        logHR += removeReassortmentEdge(edgeToRemove);        
    	
        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // HR contribution for reverse move

        if (randomlySampleAttachmentEdgeInput.get()) {
            int nPossibleSourceEdges = 0;
            for (NetworkEdge e : networkEdges) {
                if (!e.isRootEdge() && e.hasSegments.cardinality() >= 2) {
                    nPossibleSourceEdges++;
                }
            }

            logHR += Math.log(1.0/(double)nPossibleSourceEdges)
                    + Math.log(1.0/sourceEdge.getLength());
        }else {
        	double sumEdgeLengths = 0.0;
			for (NetworkEdge e : networkEdges) {
				if (!e.isRootEdge() && e.hasSegments.cardinality() >= 2) {
					sumEdgeLengths += e.getLength();
				}
			}
			logHR += Math.log(1.0 / sumEdgeLengths);
        }
 
        
        int destEdges = 0;
		for (NetworkEdge edge : networkEdges) {
			if (!edge.isRootEdge() && edge.parentNode.getHeight() > destTime
					&& edge.childNode.getHeight() <= destTime) {
				destEdges++;
			}
		}
        
        if (destEdges==0)
        	destEdges++;
        
        logHR += Math.log(1.0/destEdges);

        return logHR;
    }

    double removeReassortmentEdge(NetworkEdge edgeToRemove) {
        double logHR = 0.0;

        NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);

        // Divert segments away from chosen edge using resimulation
        BitSet segsToDivert = (BitSet) edgeToRemove.hasSegments.clone();
        logHR += divertSegmentsWithResimulation(edgeToRemoveSpouse, edgeToRemove, segsToDivert);
        
        if (divertOneSegmentInput.get()) {
            logHR += getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments, segsToDivert, 1.0/(double) edgeToRemoveSpouse.hasSegments.cardinality());
        }else {
        	logHR += getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments, segsToDivert, 0.5);
        }

        return logHR;
    }


    
    /**
     * Divert segments from sourceEdge to destEdge using resimulation approach.
     * This replaces the simple divertSegments() call with resimulation.
     */
    protected double divertSegmentsWithResimulation(NetworkEdge destEdge, NetworkEdge sourceEdge, BitSet segsToDivert) {
		double logHR = 0.0;

		BitSet segsToAdd = (BitSet) segsToDivert.clone();
		// startwith keeping a list of active edges and the segs to Add
		List<NetworkEdge> activeEdges = new ArrayList<>();
		List<BitSet> segsToAddList = new ArrayList<>();
		
		List<NetworkEdge> inactiveEdges = new ArrayList<>();
		List<BitSet> segsToRemoveList = new ArrayList<>();
		
		activeEdges.add(destEdge);
		segsToAddList.add(segsToAdd);
		
		inactiveEdges.add(sourceEdge);
		segsToRemoveList.add((BitSet) segsToDivert.clone());
		
		// Track the current time (start at the bottom of the lowest edge)
		double currentTime = destEdge.childNode.getHeight();

		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];

		getTreeNodesDown(sourceEdge, segsToDivert, treeChildNodeList);

		List<NetworkEdge> edgesAdded = new ArrayList<>();
		
		// Get network events to determine lineages at different times
		List<NetworkEvent> networkEventList = coalescentDistr.intervals.getNetworkEventList();

		logHR += divertSegmentsToAncestors(activeEdges, inactiveEdges, segsToAddList, segsToRemoveList, currentTime, edgesAdded, networkEventList, true, true);
		
		if (reconnectSegmentTrees(treeChildNodeList, destEdge, segsToDivert))
			return Double.NEGATIVE_INFINITY;
		
		cleanEmptyEdgesTopDown();
		
		return logHR;
	}
}
