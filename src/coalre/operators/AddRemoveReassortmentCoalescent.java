package coalre.operators;

import beast.base.core.Input;
import beast.base.util.Randomizer;
import coalre.distribution.CoalescentWithReassortment;
import coalre.distribution.NetworkEvent;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.statistics.NetworkStatsLogger;

import java.util.*;
import java.util.stream.Collectors;

public class AddRemoveReassortmentCoalescent extends DivertSegmentOperator {

    public Input<CoalescentWithReassortment> coalescentDistrInput = new Input<>("coalescentWithReassortment",
            "Mean of exponential used for choosing root attachment times.",
            Input.Validate.REQUIRED);
    
	public Input<Boolean> randomlySampleAttachmentEdgeInput = new Input<>("randomlySampleAttachmentEdge",
			"Randomly sample edge to attach to", true);
	
	public Input<Boolean> useDirectAncestorsOnlyInput = new Input<>("useDirectAncestorsOnly",
			"if true, only use direct ancestors of segment trees when sampling attachment edge", false);
	
	
	
    CoalescentWithReassortment coalescentDistr;

    boolean useMaxHeight = false;
    
    double addProb = 0.5;
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
        coalescentDistr = coalescentDistrInput.get();
        if (coalescentDistr.maxHeightRatioInput.get()==1.0 && coalescentDistr.redFactorInput.get()==0.0)
        	useMaxHeight = true;
    }

    @Override
    public double networkProposal() {
        double logHR;
        
        for (int i = 0; i < segmentTrees.size(); i++) {
    		segmentsChanged.set(i, false);
        }


        if (Randomizer.nextDouble()< addProb) {
        	if (useDirectAncestorsOnlyInput.get()) {
        		logHR = deltaAddRecombination();
        	}else {
        		logHR = addRecombination();
        	}
//        	System.out.println(network +"\n");
//        	System.exit(0);
        	logHR += Math.log((1.0 - addProb)/addProb);
        }else {
        	try {
				if (useDirectAncestorsOnlyInput.get())
					logHR = deltaRemoveRecombination();
				else
					logHR = removeRecombination();
                logHR += Math.log(addProb/(1.0 - addProb));
//                System.out.println("remove: " + logHR);
        	} catch (Exception e) {
        		return Double.NEGATIVE_INFINITY;
        	}
        }
        
        return logHR;
    }

    double addRecombination() {
        double logHR = 0.0;

        List<Integer> possibleSourceEdges = new ArrayList<>();
        if (useMaxHeight) {
    		double lociMRCA = NetworkStatsLogger.getSecondHighestLociMRCA(networkInput.get());

        	// If we use max height, we can only add reassortment edges to the root edge
			for (int i = 0; i < networkEdges.size(); i++) {
				NetworkEdge edge = networkEdges.get(i);
				if (!edge.isRootEdge() 
						&& edge.hasSegments.cardinality() >= 2
						&& edge.childNode.getHeight() < lociMRCA) {
					possibleSourceEdges.add(i);
				}
			}
        }else {
			for (int i = 0; i < networkEdges.size(); i++) {
				NetworkEdge edge = networkEdges.get(i);
				if (!edge.isRootEdge() && edge.hasSegments.cardinality() >= 2) {
					possibleSourceEdges.add(i);
				}
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
					
//					possibleSourceEdges.stream().mapToDouble(e -> e.getLength()).sum();			
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
        // keep only those that coexist at the time of maxHeight
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
        
        int nRemovableEdges;
//		if (divertOneSegmentInput.get()) {
//			nRemovableEdges= (int) networkEdges.stream().filter(e -> !e.isRootEdge())
//				.filter(e -> e.hasSegments.cardinality() == 1).filter(e -> e.childNode.isReassortment())
//				.filter(e -> e.parentNode.isCoalescence()).count();
//		} else {
			nRemovableEdges= (int) networkEdges.stream().filter(e -> !e.isRootEdge())
					.filter(e -> e.hasSegments.cardinality() >= 1).filter(e -> e.childNode.isReassortment())
					.filter(e -> e.parentNode.isCoalescence()).count();

//		}


//        // HR contribution for reverse move
//        int nRemovableEdges = (int) networkEdges.stream()
//                .filter(e -> !e.isRootEdge())
//                .filter(e -> e.hasSegments.cardinality()>=1)
//                .filter(e -> e.childNode.isReassortment())
//                .filter(e -> e.parentNode.isCoalescence())
//                .count();
        
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
        if (divertOneSegmentInput.get()) {
        	double prob = 1./sourceEdge.hasSegments.cardinality();
            BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments, prob);
            logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments, segsToDivert, prob);
            logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
        }else {
	        BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
            logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments);
            logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
        }
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
        
//        List<NetworkEdge> removableEdges = networkEdges.stream()
//                .filter(e -> !e.isRootEdge())
//                .filter(e -> e.hasSegments.cardinality()>=1)
//                .filter(e -> e.childNode.isReassortment())
//                .filter(e -> e.parentNode.isCoalescence())
//                .collect(Collectors.toList());
        
        
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
    	
//    	double logHRBefore = logHR;
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
//    	System.out.println(logHR-logHRBefore);
    	
        // Remove reassortment edge
        logHR += removeReassortmentEdge(edgeToRemove);        
    	
        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // HR contribution for reverse move


        if (randomlySampleAttachmentEdgeInput.get()) {
            int nPossibleSourceEdges = 0;
            if (useMaxHeight) {
            	double lociMRCA = NetworkStatsLogger.getSecondHighestLociMRCA(networkInput.get());
				for (NetworkEdge e : networkEdges) {
					if (!e.isRootEdge() && e.hasSegments.cardinality() >= 2 && e.childNode.getHeight() < lociMRCA) {
						nPossibleSourceEdges++;
					}
				}
            }else {
	            for (NetworkEdge e : networkEdges) {
	                if (!e.isRootEdge() && e.hasSegments.cardinality() >= 2) {
	                    nPossibleSourceEdges++;
	                }
	            }
            }

            logHR += Math.log(1.0/(double)nPossibleSourceEdges)
                    + Math.log(1.0/sourceEdge.getLength());
        }else {
        	double sumEdgeLengths = 0.0;
        	if (useMaxHeight) {
        		double lociMRCA = NetworkStatsLogger.getSecondHighestLociMRCA(networkInput.get());
        		for (NetworkEdge e : networkEdges) {
        			if (!e.isRootEdge() && e.hasSegments.cardinality() >= 2 && e.childNode.getHeight() < lociMRCA) {
        				sumEdgeLengths += e.getLength();
        			}
        		}
        	}else {
				for (NetworkEdge e : networkEdges) {
					if (!e.isRootEdge() && e.hasSegments.cardinality() >= 2) {
						sumEdgeLengths += e.getLength();
					}
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


        NetworkNode nodeToRemove = edgeToRemove.childNode;
        NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
        NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;
		networkEdges.remove(edgeToRemove);
		networkEdges.remove(edgeToRemoveSpouse);

        // Divert segments away from chosen edge
        BitSet segsToDivert = (BitSet) edgeToRemove.hasSegments.clone();
        logHR +=  divertSegments(edgeToRemoveSpouse, edgeToRemove, segsToDivert);
        if (divertOneSegmentInput.get()) {
//        	System.out.println(getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments, segsToDivert, 1.0/(double) edgeToRemoveSpouse.hasSegments.cardinality()));
            logHR += getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments, segsToDivert, 1.0/(double) edgeToRemoveSpouse.hasSegments.cardinality());
        }else {
        	logHR += getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments);
        }

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

        networkEdges.remove(secondNodeToRemove.getParentEdges().get(0));
        if (secondNodeToRemove.getParentEdges().get(0).isRootEdge()) {
            network.setRootEdge(secondEdgeToExtend);

        } else {
            NetworkEdge secondNodeToRemoveParentEdge = secondNodeToRemove.getParentEdges().get(0);
            NetworkNode secondNodeToRemoveParent = secondNodeToRemoveParentEdge.parentNode;
            secondNodeToRemoveParent.removeChildEdge(secondNodeToRemoveParentEdge);
            secondNodeToRemove.removeParentEdge(secondNodeToRemoveParentEdge);

            secondNodeToRemoveParent.addChildEdge(secondEdgeToExtend);
            networkEdges.remove(secondNodeToRemoveParentEdge);
        }

        return logHR;
    }

    
    // Delta-based add recombination move
    double deltaAddRecombination() {
        double logHR = 0.0;

        // 1. Select source edge and time (same as original)
        List<Integer> possibleSourceEdges = new ArrayList<>();
        if (useMaxHeight) {
            double lociMRCA = NetworkStatsLogger.getSecondHighestLociMRCA(networkInput.get());
            for (int i = 0; i < networkEdges.size(); i++) {
                NetworkEdge edge = networkEdges.get(i);
                if (!edge.isRootEdge() 
                        && edge.hasSegments.cardinality() >= 2
                        && edge.childNode.getHeight() < lociMRCA) {
                    possibleSourceEdges.add(i);
                }
            }
        } else {
            for (int i = 0; i < networkEdges.size(); i++) {
                NetworkEdge edge = networkEdges.get(i);
                if (!edge.isRootEdge() && edge.hasSegments.cardinality() >= 2) {
                    possibleSourceEdges.add(i);
                }
            }
        }
        
        NetworkEdge sourceEdge;
        double sourceTime;
        
        if (randomlySampleAttachmentEdgeInput.get()) {
            sourceEdge = networkEdges.get(possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size())));
            sourceTime = Randomizer.nextDouble() * sourceEdge.getLength() + sourceEdge.childNode.getHeight();
            logHR -= Math.log(1.0/(double)possibleSourceEdges.size())
                    + Math.log(1.0/sourceEdge.getLength());
        } else {
            double sumEdgeLengths = 0.0;
            for (Integer i : possibleSourceEdges) {
                sumEdgeLengths += networkEdges.get(i).getLength();
            }
            double randomEdgeLength = Randomizer.nextDouble() * sumEdgeLengths;
            logHR -= Math.log(1.0/sumEdgeLengths);
            double passedLength = 0;
            sourceEdge = null;
            sourceTime = -1.0;
            for (Integer i : possibleSourceEdges) {
                passedLength += networkEdges.get(i).getLength();
                if (passedLength > randomEdgeLength) {
                    sourceEdge = networkEdges.get(i);
                    sourceTime = passedLength-randomEdgeLength + sourceEdge.childNode.getHeight();
                    break;
                }
            }
        }

        // 2. Sample attachment time using coalescent intervals (same as original)
        List<NetworkEvent> networkEventList = coalescentDistr.intervals.getNetworkEventList();
        double currTime = sourceTime;
        NetworkEvent prevEvent = null;
        double attachmentTime = 0.0;
        
        for (NetworkEvent event : networkEventList) {
            if (event.time > currTime) {
                double rate = 0.1 * prevEvent.lineages;
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
        
        if (attachmentTime > network.getRootEdge().childNode.getHeight()) {
            double currentTransformedTime = coalescentDistr.populationFunction.getIntensity(network.getRootEdge().childNode.getHeight());
            double transformedTimeToNextCoal = Randomizer.nextExponential(0.5);
            attachmentTime = coalescentDistr.populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime);
            
            logHR -= -0.1 * coalescentDistr.populationFunction.getIntegral(network.getRootEdge().childNode.getHeight(), attachmentTime);
            logHR -= Math.log(0.1/coalescentDistr.populationFunction.getPopSize(attachmentTime));
        }

        double destTime = attachmentTime;

        List<NetworkEdge> targetEdges = new ArrayList<>();
        List<NetworkEdge> visitedEdges = new ArrayList<>();
        getTargetEdges(sourceEdge, destTime, targetEdges, visitedEdges);
        
        // only keep unique target Edges
        targetEdges = targetEdges.stream()
                .distinct()
                .collect(Collectors.toList());


        
        if (targetEdges.isEmpty()) {
            return Double.NEGATIVE_INFINITY;
        }

        NetworkEdge destEdge = targetEdges.get(Randomizer.nextInt(targetEdges.size()));
        logHR -= Math.log(1.0/targetEdges.size());

        if (!destEdge.isRootEdge() && destEdge.parentNode.getHeight() < sourceTime)
            return Double.NEGATIVE_INFINITY;

        // 4. Create new reassortment edge
        logHR += addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);
        
        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        
        
        // 1. Select removable edge - only those that could have been added by deltaAddRecombination
        List<Integer> removableEdges = new ArrayList<>();
        for (int i = 0; i < networkEdges.size(); i++) {
            NetworkEdge edge = networkEdges.get(i);
            if (!edge.isRootEdge() && edge.childNode.isReassortment()
                    && edge.parentNode.isCoalescence() && edge.hasSegments.cardinality() >= 1) {
                
                NetworkEdge potentialSourceEdge = getSpouseEdge(edge);
                double trueDestTime = edge.parentNode.getHeight();
                NetworkEdge trueDestEdge = getSisterEdge(edge);
                
                // Check if destEdge is reachable from sourceEdge via delta traversal
                List<NetworkEdge> reachableEdges = new ArrayList<>();
                visitedEdges = new ArrayList<>();
                getTargetEdges(potentialSourceEdge, trueDestTime, reachableEdges, visitedEdges);
                
//                System.out.println("Checking edge: " + trueDestTime + " " + trueDestEdge.childNode.getHeight() + " - " + trueDestEdge.parentNode.getHeight());
//                
//                for (NetworkEdge e : reachableEdges)
//                	System.out.println(e.childNode.getHeight() + " - " + e.parentNode.getHeight());
//                
                if (reachableEdges.contains(trueDestEdge)) {
                    removableEdges.add(i);
                }
            }
        }
        
        // only keep unique edges
        Collections.sort(removableEdges);
        removableEdges = removableEdges.stream().distinct().collect(Collectors.toList());
        logHR += Math.log(1.0/removableEdges.size());
        return logHR;
    }

    // Delta-based remove recombination move
    double deltaRemoveRecombination() {
        double logHR = 0.0;
        
        // 1. Select removable edge - only those that could have been added by deltaAddRecombination
        List<Integer> removableEdges = new ArrayList<>();
        for (int i = 0; i < networkEdges.size(); i++) {
            NetworkEdge edge = networkEdges.get(i);
            if (!edge.isRootEdge() && edge.childNode.isReassortment()
                    && edge.parentNode.isCoalescence() && edge.hasSegments.cardinality() >= 1) {
                
                NetworkEdge potentialSourceEdge = getSpouseEdge(edge);
                double trueDestTime = edge.parentNode.getHeight();
                NetworkEdge trueDestEdge = getSisterEdge(edge);
                
                
                // Check if destEdge is reachable from sourceEdge via delta traversal
                List<NetworkEdge> reachableEdges = new ArrayList<>();
                List<NetworkEdge> visitedEdges = new ArrayList<>();
                getTargetEdges(potentialSourceEdge, trueDestTime, reachableEdges, visitedEdges);
                
//                System.out.println("Checking edge: " + trueDestTime + " " + trueDestEdge.childNode.getHeight() + " - " + trueDestEdge.parentNode.getHeight());
//                
//                for (NetworkEdge e : reachableEdges)
//                	System.out.println(e.childNode.getHeight() + " - " + e.parentNode.getHeight());
//                
                if (reachableEdges.contains(trueDestEdge) || reachableEdges.contains(trueDestEdge.parentNode.getParentEdges().get(0))) {
                    removableEdges.add(i);
                }
            }
        }
        // only keep unique edges
        Collections.sort(removableEdges);
        removableEdges = removableEdges.stream().distinct().collect(Collectors.toList());
        
        
        if (removableEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        NetworkEdge edgeToRemove = networkEdges.get(removableEdges.get(Randomizer.nextInt(removableEdges.size())));
        logHR -= Math.log(1.0/(removableEdges.size()));

        // 2. Extract source/dest information
        double sourceTime = edgeToRemove.childNode.getHeight();
        NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
        NetworkEdge destEdge = getSisterEdge(edgeToRemove);
        if (destEdge.childNode == edgeToRemove.childNode)
            destEdge = sourceEdge;
        double destTime = edgeToRemove.parentNode.getHeight();

        // 3. Time density contribution (same as original)
        List<NetworkEvent> networkEventList = coalescentDistr.intervals.getNetworkEventList();
        double currTime = sourceTime;
        NetworkEvent prevEvent = null;
        
        for (NetworkEvent event : networkEventList) {
            if (event.time > currTime) {
                double rate = 0.5 * (prevEvent.lineages - 1);
                if (destTime <= event.time) {
                    logHR += -rate * coalescentDistr.populationFunction.getIntegral(currTime, destTime);
                    logHR += Math.log(rate / coalescentDistr.populationFunction.getPopSize(destTime));
                    break;
                }
                logHR += -rate * coalescentDistr.populationFunction.getIntegral(currTime, event.time);
                currTime = event.time;
            }
            prevEvent = event;
        }

        // 4. Remove reassortment edge
        logHR += removeReassortmentEdge(edgeToRemove);
        
        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // 5. Reverse move probability components
        // 5a. Source edge/time selection
        if (randomlySampleAttachmentEdgeInput.get()) {
            int nPossibleSourceEdges = 0;
            if (useMaxHeight) {
                double lociMRCA = NetworkStatsLogger.getSecondHighestLociMRCA(networkInput.get());
                for (NetworkEdge e : networkEdges) {
                    if (!e.isRootEdge() && e.hasSegments.cardinality() >= 2 && e.childNode.getHeight() < lociMRCA) {
                        nPossibleSourceEdges++;
                    }
                }
            } else {
                for (NetworkEdge e : networkEdges) {
                    if (!e.isRootEdge() && e.hasSegments.cardinality() >= 2) {
                        nPossibleSourceEdges++;
                    }
                }
            }
            logHR += Math.log(1.0/(double)nPossibleSourceEdges) + Math.log(1.0/sourceEdge.getLength());
        } else {
            double sumEdgeLengths = 0.0;
            if (useMaxHeight) {
                double lociMRCA = NetworkStatsLogger.getSecondHighestLociMRCA(networkInput.get());
                for (NetworkEdge e : networkEdges) {
                    if (!e.isRootEdge() && e.hasSegments.cardinality() >= 2 && e.childNode.getHeight() < lociMRCA) {
                        sumEdgeLengths += e.getLength();
                    }
                }
            } else {
                for (NetworkEdge e : networkEdges) {
                    if (!e.isRootEdge() && e.hasSegments.cardinality() >= 2) {
                        sumEdgeLengths += e.getLength();
                    }
                }
            }
            logHR += Math.log(1.0 / sumEdgeLengths);
        }


        List<NetworkEdge> reverseTargetEdges = new ArrayList<>();
        List<NetworkEdge> visitedEdges = new ArrayList<>();
        getTargetEdges(sourceEdge, destTime, reverseTargetEdges, visitedEdges);
        reverseTargetEdges = reverseTargetEdges.stream()
                .distinct()
                .collect(Collectors.toList());

        
        int numCandidates = Math.max(1, reverseTargetEdges.size());
        logHR += Math.log(1.0 / (double) numCandidates);

        return logHR;
    }

    
    // Traverse up the network (forward in time) - only direct ancestors
    private void getTargetEdges(NetworkEdge currentEdge, double destTime, 
                                List<NetworkEdge> targetEdges, List<NetworkEdge> visitedEdges) {
        
        // Check if this edge is already in target edges - prevent exponential growth
        if (visitedEdges.contains(currentEdge)) {
            return;
        }
        visitedEdges.add(currentEdge);
        
        if (currentEdge.isRootEdge()) {
            // Reached root
            targetEdges.add(currentEdge);
            return;
        }
        
		if (currentEdge.childNode.getHeight() <= destTime && currentEdge.parentNode.getHeight() >= destTime) {
			// Current edge spans destTime - add to target edges
			targetEdges.add(currentEdge);
            // if the parentNode is a coalescent event, check the sister edge as well
//            if (currentEdge.parentNode.isCoalescence()) {
//                NetworkEdge sisterEdge = getSisterEdge(currentEdge);
//                if (sisterEdge.childNode.getHeight() <= destTime && sisterEdge.parentNode.getHeight() >= destTime) {
//                    targetEdges.add(sisterEdge);
//                }
//            }
			return;
		}else {
            for (NetworkEdge e : currentEdge.parentNode.getParentEdges()) {
            	getTargetEdges(e, destTime, targetEdges, visitedEdges);
            }

		}
    }
    
 }
