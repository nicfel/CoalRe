package coalre.operators;

import beast.base.core.Input;
import beast.base.util.Randomizer;
import coalre.distribution.CoalescentWithReassortment;
import coalre.distribution.NetworkEvent;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.*;
import java.util.stream.Collectors;

public class AddRemoveReassortmentCoalescent extends DivertSegmentOperator {

    public Input<CoalescentWithReassortment> coalescentDistrInput = new Input<>("coalescentWithReassortment",
            "Mean of exponential used for choosing root attachment times.",
            Input.Validate.REQUIRED);
    
	public Input<Boolean> randomlySampleAttachmentEdgeInput = new Input<>("randomlySampleAttachmentEdge",
			"Randomly sample edge to attach to", true);
	

    CoalescentWithReassortment coalescentDistr;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        coalescentDistr = coalescentDistrInput.get();
    }

    @Override
    public double networkProposal() {
        double logHR;
        
//        System.out.println(network);
        for (int i = 0; i < segmentTrees.size(); i++) {
    		segmentsChanged.set(i, false);
        }


        if (Randomizer.nextBoolean()) {
            logHR = addRecombination();
        }else {
        	try {
        		logHR = removeRecombination();
        	} catch (Exception e) {
        		return Double.NEGATIVE_INFINITY;
        	}
        }
//        System.out.println(network);
        return logHR;
    }

    double addRecombination() {
        double logHR = 0.0;
//        System.out.println(network);

        List<Integer> possibleSourceEdges = new ArrayList<>();
		for (int i = 0; i < networkEdges.size(); i++) {
			NetworkEdge edge = networkEdges.get(i);
			if (!edge.isRootEdge() && edge.hasSegments.cardinality() >= 2) {
				possibleSourceEdges.add(i);
			}
		}
        
//        List<NetworkEdge> possibleSourceEdges = networkEdges.stream()
//                .filter(e -> !e.isRootEdge())
//                .filter(e -> e.hasSegments.cardinality()>=2)
//                .collect(Collectors.toList());
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

        // HR contribution for reverse move
        int nRemovableEdges = (int) networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=1)
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.parentNode.isCoalescence())
                .count();
        
        logHR += Math.log(1.0/nRemovableEdges);
//        System.out.println(network);
//        System.exit(0);
        return logHR;
    }
    
    double addReassortmentEdge(NetworkEdge sourceEdge, double sourceTime,
                               NetworkEdge destEdge, double destTime) {

        double logHR = 0.0;


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
        if (divertOneSegmentInput.get()) {
            BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments, 0.1);
            
            for (int i = 0; i < segmentTrees.size(); i++) {
            	if (segsToDivert.get(i))
            		segmentsChanged.set(i, true);            	
            }

            
            logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments, segsToDivert, 0.1);
            logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);

//            logHR -= addSegmentsToAncestors(reassortmentEdge, segsToDivert);
//            logHR += removeSegmentsFromAncestors(newEdge1, segsToDivert);
            // get the logHR in the other direction
//            logHR += getLogConditionedSubsetProb(sourceEdge.hasSegments, segsToDivert, 0.1);
        }else {
	        BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
            for (int i = 0; i < segmentTrees.size(); i++) {
            	if (segsToDivert.get(i))
            		segmentsChanged.set(i, true);            	
            }
            logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments);
            logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);

//	        logHR -= addSegmentsToAncestors(reassortmentEdge, segsToDivert);
//	        logHR += removeSegmentsFromAncestors(newEdge1, segsToDivert);
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
//            		(int)networkEdges.stream()
//                    .filter(e -> !e.isRootEdge())
//                    .filter(e -> e.hasSegments.cardinality()>=2)
//                    .count();

            logHR += Math.log(1.0/(double)nPossibleSourceEdges)
                    + Math.log(1.0/sourceEdge.getLength());
        }else {
        	double sumEdgeLengths = 0.0;
			for (NetworkEdge e : networkEdges) {
				if (!e.isRootEdge() && e.hasSegments.cardinality() >= 2) {
					sumEdgeLengths += e.getLength();
				}
			}
//			List<NetworkEdge> possibleSourceEdges = networkEdges.stream().filter(e -> !e.isRootEdge())
//					.filter(e -> e.hasSegments.cardinality() >= 2).collect(Collectors.toList());
//			double sumEdgeLengths = possibleSourceEdges.stream().mapToDouble(e -> e.getLength()).sum();
			logHR += Math.log(1.0 / sumEdgeLengths);
        }
 
        
        int destEdges = 0;
		for (NetworkEdge edge : networkEdges) {
			if (!edge.isRootEdge() && edge.parentNode.getHeight() > destTime
					&& edge.childNode.getHeight() <= destTime) {
				destEdges++;
			}
		}
//        // keep only those that coexist at the time of maxHeight
//        List<NetworkEdge> destEdges = networkEdges.stream()
//                .filter(e -> !e.isRootEdge())
//                .filter(e -> e.parentNode.getHeight()>destTime)
//                .filter(e -> e.childNode.getHeight()<=destTime)
//               .collect(Collectors.toList());
        
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

        // Divert segments away from chosen edge
        BitSet segsToDivert = (BitSet) edgeToRemove.hasSegments.clone();
        for (int i = 0; i < segmentTrees.size(); i++) {
        	if (segsToDivert.get(i))
        		segmentsChanged.set(i, true);            	
        }

        logHR +=  divertSegments(edgeToRemoveSpouse, edgeToRemove, segsToDivert);

        
        if (divertOneSegmentInput.get()) {
            logHR += getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments, segsToDivert, 0.1);
		} else {
			logHR += getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments);
		}
        
        networkEdges.remove(edgeToRemove);
        networkEdges.remove(edgeToRemoveSpouse);


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
        }
//        if (!networkTerminatesAtMRCA())
//            return Double.NEGATIVE_INFINITY;

        return logHR;
    }
    
    
    private BitSet getRandomConditionedSubset(BitSet sourceSegments, double binomProb) {

        if (sourceSegments.cardinality() < 2) {
            return null;
        }

        BitSet destSegments = new BitSet();

        do {
            destSegments.clear();

            for (int segIdx = sourceSegments.nextSetBit(0);
                 segIdx != -1;
                 segIdx = sourceSegments.nextSetBit(segIdx + 1)) {

                // Now pick each segment with probability 0.1
                if (Randomizer.nextDouble() < binomProb) {
                    destSegments.set(segIdx);
                }
            }

        } while (destSegments.cardinality() == 0
                || destSegments.cardinality() == sourceSegments.cardinality());
//        System.out.println(destSegments + " " + sourceSegments);
        return destSegments;
    }
    
    
    private double getLogConditionedSubsetProb(BitSet sourceSegments, BitSet chosenSubset, double binomProb) {

        if (sourceSegments.cardinality() < 2) {
            return Double.NEGATIVE_INFINITY;
        }

        // Let n be the total number of possible segments in sourceSegments,
        // and k be the number actually chosen in the subset.
        int n = sourceSegments.cardinality();
        int k = chosenSubset.cardinality();


        // Fully correct conditional log-probability:
        //    log( p^k * (1-p)^(n-k) )  -  log( 1 - p^n - (1-p)^n )
        // =  k*log(p) + (n-k)*log(1-p) - log( 1 - p^n - (1-p)^n )
        double logProb = k * Math.log(binomProb)
                       + (n - k) * Math.log(1.0 - binomProb)
                       - Math.log(1.0 - Math.pow(binomProb, n) - Math.pow(1.0 - binomProb, n));

        return logProb;
    }

    
    
//    /**
//     * automatic parameter tuning *
//     */
//    @Override
//    public void optimize(final double logAlpha) {
//        if (optimiseInput.get()) {
//            double delta = calcDelta(logAlpha);
//            delta += Math.log(alpha);
//            final double f = Math.exp(delta);
//            if( alpha > 0 ) {
//                final RecombinationNetwork network = networkInput.get();
//                final double h = network.getRootEdge().childNode.getHeight();
//                final double k = Math.log(network.getLeafNodes().size()) / Math.log(2);
//                final double lim = (h / k) * alpha;
//                if( f <= lim ) {
//                	alpha = f;
//                }
//            } else {
//            	alpha = f;
//            }
//        }
//    }

 }
