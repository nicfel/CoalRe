package coalre.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import coalre.distribution.CoalescentWithReassortment;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.*;
import java.util.stream.Collectors;

public class AddRemoveReassortmentExponential extends AddRemoveReassortment {

	public Input<Boolean> randomlySampleAttachmentEdgeInput = new Input<>("randomlySampleAttachmentEdge",
			"Randomly sample edge to attach to", true);


    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<Boolean> gaussianInput = new Input<>("gaussian", "Gaussian (=true=default) or uniform delta", true);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, " +
            "i.e. -1. (when positive, gets multiplied by network-height/log2(n-taxa).", -1.0);

	
	CoalescentWithReassortment coalescentDistr;

	
    // shadows size
    double size;
    private double limit;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
        size = sizeInput.get();
        limit = limitInput.get();
        
		if (segmentTrees.isEmpty()) {
			throw new IllegalStateException("No segment trees available for recombination.");
		}

	}

	@Override
	public double networkProposal() {
		double logHR=0.0;
		if (Randomizer.nextBoolean()) {
			logHR = addRecombination();
		} else {
//			try {
				logHR = removeRecombination();
//				System.out.println("remove " + logHR);
//			} catch (Exception e) {
//				return Double.NEGATIVE_INFINITY;
//			}
		}		
		return logHR;
	}

	double addRecombination() {			
		double logHR = 0.0;
		// pick a random edge
		List<NetworkEdge> possibleSourceEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() >= 2)
				.collect(Collectors.toList());

		NetworkEdge sourceEdge = possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size()));
		double sourceTime = Randomizer.nextDouble() * sourceEdge.getLength() + sourceEdge.childNode.getHeight();

		logHR -= Math.log(1.0 / (double) possibleSourceEdges.size()) + Math.log(1.0 / sourceEdge.getLength());

		// pick a random segment on the source edge
		int segment = Randomizer.nextInt(network.getSegmentCount());
		while (!sourceEdge.hasSegments.get(segment)) {
			segment = Randomizer.nextInt(network.getSegmentCount());
		}
		logHR -= Math.log(1.0 / sourceEdge.hasSegments.cardinality()); // pick random segment
		
		// follow the edge down until we find the next node with the segment
		NetworkEdge sourceStartEdge = getSegmentNodeIndex(segment, sourceEdge);
				
		double difference = sourceStartEdge.childNode.getHeight() - sourceTime;
        double delta = Math.abs(getDelta(difference));        
        logHR -= getLogDeltaProb(delta, difference); //pick delta

        // get all potential reattachment Edges
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();
		// calculate the max height of the new child
        
		if (sourceTime + delta > sourceEdge.parentNode.getHeight()) {
			// skip source edge
			delta -= (sourceEdge.parentNode.getHeight()-sourceTime);
			for (NetworkEdge e : sourceEdge.parentNode.getParentEdges()) {
				if (e.hasSegments.get(segment)) {
			    	targetEdges.putAll(getTargetEdgesUp(e, delta, sourceTime, segment));
				}
			}
		} else {
	    	targetEdges.put(sourceEdge, sourceTime + delta); // reattach to source edge)
		}

    	
    	logHR -= Math.log(1.0/targetEdges.size());

        List<NetworkEdge> potentialNewTargets = targetEdges.keySet().stream().collect(Collectors.toList());        
        NetworkEdge destEdge = potentialNewTargets.get(Randomizer.nextInt(potentialNewTargets.size()));
        double destTime = targetEdges.get(destEdge);
        
        		
		// Create new reassortment edge
		logHR += addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime, segment);	
		
		if (logHR == Double.NEGATIVE_INFINITY)
			return Double.NEGATIVE_INFINITY;

		int nRemovableEdges = (int) networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.childNode.isReassortment())
				.filter(e -> e.parentNode.isCoalescence())
				.count();
		
		// pick the newly added edge to remove
		logHR += Math.log(1.0/nRemovableEdges);
//	    BitSet coalSegments = (BitSet) destEdge.hasSegments.clone();
//	    coalSegments.and(getSisterEdge(destEdge).hasSegments);
////    	logHR += Math.log(1.0 / coalSegments.cardinality()); // pick the same segment in the reverse move
		return logHR;
	}



	double removeRecombination() {
	    double logHR = 0.0;

	    List<Integer> removableEdges = new ArrayList<>();
	    for (int i = 0; i < networkEdges.size(); i++) {
	        NetworkEdge edge = networkEdges.get(i);
	        if (!edge.isRootEdge() 
	        		&& edge.childNode.isReassortment() 
	        		&& edge.parentNode.isCoalescence()) {
	            removableEdges.add(i);
	        }
	    }
	    
	    if (removableEdges.isEmpty()) {
	        return Double.NEGATIVE_INFINITY;
	    }
	    
	    NetworkEdge edgeToRemove = networkEdges.get(removableEdges.get(Randomizer.nextInt(removableEdges.size())));
	    logHR -= Math.log(1.0 / removableEdges.size());
	    
	    
	    // pick a random segment to later compute the distance delta and reattachment times from
	    int segment =-1;
	    while (segment==-1) {
	    	segment = Randomizer.nextInt(network.getSegmentCount()); 
	    	if (!edgeToRemove.hasSegments.get(segment))
	    		segment = -1; // segment not available, try again	    	   
	    }
    	logHR -= Math.log(1.0 / edgeToRemove.hasSegments.cardinality()); // pick random segment

	    
	    // keep track of the timings
	    double destTime = edgeToRemove.parentNode.getHeight();	    
	    double sourceTime = edgeToRemove.childNode.getHeight();
	    NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
	    NetworkEdge destEdge = getSisterEdge(edgeToRemove);
	    
	    if (destEdge.childNode == edgeToRemove.childNode)
	        destEdge = sourceEdge;
	    
	    
//	    System.out.println(network);
//	    double prevRootHeight = segmentTree.getRoot().getHeight();
//	    boolean wasRoot = startNode.getParent().isRoot();
	    logHR += removeReassortmentEdge(edgeToRemove, segment);
//	    System.out.println(network);
	    if (logHR == Double.NEGATIVE_INFINITY) {
	        return Double.NEGATIVE_INFINITY;
	    }
	    
	    // calculate the distance from source time to dest time
	    
	    
	    
	    
	    
//
//	    double delta = 0.0;
//	    // if start node or destNodebefore are root, do something else
//	    if (wasRoot) {
//	    	delta = startNode.getParent().getHeight() - prevRootHeight;
//	    } else if (startNode.getParent().isRoot()) {
//			// is newly created root
//	    	delta = oldDestTime - startNode.getParent().getHeight();
//		}else {
//		    double mrcaTime = getMRCAtime(startNode, destNodeBefore);
//		    double pathLength = 2 * mrcaTime - startNode.getHeight() - oldDestTime;
//		    delta = pathLength - startNode.getLength();
//		}
//	    
//	    // compute the implied delta
//	    delta = Math.abs(delta);
//	    
//		
//		Node otherChild = startNode.getParent().getLeft() == startNode ? startNode.getParent().getRight() : startNode.getParent().getLeft();
//		double minHeight = startNode.getHeight();
//
//		Map<Node, Double> targetEdges = new HashMap<>();
//    	targetEdges.putAll(getTargetEdgesUp(startNode.getParent(), delta, minHeight, segment));
//
//    	// pick a random segment
//    	logHR += Math.log(1.0 / network.getSegmentCount());
//    	// pick any node
//    	logHR += Math.log(1.0 / (segmentTree.getNodeCount() - 1));		
//    	// calculate probability of the required delta
//		logHR += getLogDeltaProb(delta, startNode.getLength());	
//		// pick a target edge
//    	logHR += Math.log(1.0 / targetEdges.size());
//    	
//    	// contribution for the source time
//        double timeDiff = Math.min(startNode.getParent().getHeight(), oldDestTime) - startNode.getHeight();
//        logHR += Math.log(1.0 / timeDiff);

	    return logHR;
	}

	private double getMRCAtime(Node n1, Node n2) {
		if (n1==n2)
			return n1.getHeight();
		if (n1.getHeight() > n2.getHeight()) {
			return getMRCAtime(n1, n2.getParent());
		}else {
			return getMRCAtime(n1.getParent(), n2);
		}
	}

	
	double addReassortmentEdge(NetworkEdge sourceEdge, double sourceTime, NetworkEdge destEdge, double destTime,
			int segment) {

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

		BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments, segment);
//		System.out.println(segsToDivert);
		logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments);		
		logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
		

		networkEdges.add(reassortmentEdge);
		networkEdges.add(newEdge1);
		networkEdges.add(newEdge2);
		return logHR;
	}


	
	
	double removeReassortmentEdge(NetworkEdge edgeToRemove, int segment) {
		double logHR = 0.0;

		NetworkNode nodeToRemove = edgeToRemove.childNode;
		NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
		NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;
		networkEdges.remove(edgeToRemove);
		networkEdges.remove(edgeToRemoveSpouse);

		// Divert segments away from chosen edge
		BitSet segsToDivert = (BitSet) edgeToRemove.hasSegments.clone();

		logHR += divertSegments(edgeToRemoveSpouse, edgeToRemove, segsToDivert);
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

    
    private Map<NetworkEdge, Double> getTargetEdgesUp(NetworkEdge edge, double delta, double minHeight, int segment) {
    	
    	// initalize map
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();
        if (!edge.hasSegments.get(segment))
        	return targetEdges; // no valid targets
        
        
  		if (edge.isRootEdge()) {
			double newHeight = edge.childNode.getHeight() + delta;
			targetEdges.put(edge, newHeight);			
			return targetEdges;
		}        
  			
        // get the difference between the current height and the new height
        delta -= edge.getLength();

		if (delta <= 0.0) {
			double newHeight = edge.childNode.getHeight() + delta + edge.getLength();

			targetEdges.put(edge, newHeight);			
			return targetEdges; 
		} else {
			// proceed to the next event
			if (edge.parentNode.isReassortment()) {
				// if reassortment, follow both parents
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(0), delta, minHeight, segment));
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(1), delta, minHeight, segment));				
			} else {
				// else it is a coalescent event
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(0), delta, minHeight, segment));
				targetEdges.putAll(getTargetEdgesDown(getSisterEdge(edge), delta, minHeight, segment));
			}
		}
		return targetEdges;
	}

	
    private Map<NetworkEdge, Double> getTargetEdgesDown(NetworkEdge edge, double delta, double minHeight, int segment) {
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();
        if (!edge.hasSegments.get(segment))
        	return targetEdges; // no valid targets

        // get the difference between the current height and the new height
        delta -= edge.getLength();
//        if (visitedNodes.contains(edge)) {
//        	return targetEdges; // do not revisit nodes
//        }
//        visitedNodes.add(edge);
        if (delta <= 0.0) {
        	double newHeight = edge.parentNode.getHeight() - delta - edge.getLength();
            if (newHeight < minHeight)
            	return targetEdges; // no valid targets
            
			targetEdges.put(edge, edge.parentNode.getHeight() - (delta + edge.getLength()));
			return targetEdges; 
		} else {
			if (edge.childNode.getHeight() < minHeight) {
				return targetEdges; // no valid targets
			}

			// proceed to the next event
			if (edge.childNode.isLeaf()) {
				return targetEdges; // no valid targets
			}else if (edge.childNode.isReassortment()) {
				// if reassortment, follow the other parent up and the child down
//				targetEdges.putAll(getTargetEdgesUp(getSpouseEdge(edge), delta, minHeight, iEdge, visitedNodes));
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(0), delta, minHeight, segment));
			} else {
				// if coalescent event, follow both children down
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(0), delta, minHeight, segment));
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(1), delta, minHeight, segment));			
			}
		}
        return targetEdges;
	}

	
    private boolean actualCoalescent(NetworkEdge e) {
		int segment = e.hasSegments.nextSetBit(0);
		return getSisterEdge(e).hasSegments.get(segment);
	}

	private NetworkEdge getParentEdge(NetworkEdge edge, double time, int segment) {
		if (edge.isRootEdge())
			return edge;
		
		if (edge.parentNode.getHeight() > time) {
            return edge;
        }else if (edge.parentNode.isCoalescence()) {
			return getParentEdge(edge.parentNode.getParentEdges().get(0), time, segment);
		}else {
			if (edge.parentNode.getParentEdges().get(0).hasSegments.get(segment)) {
				return getParentEdge(edge.parentNode.getParentEdges().get(0), time, segment);
			}else {
				return getParentEdge(edge.parentNode.getParentEdges().get(1), time, segment);
			}
		}
				
	}

	private Optional<NetworkEdge> findNetworkEdge(int sourceNodeIndex, int segment, double height) {		
		return networkEdges.stream()
				.filter(e -> e.childNode.getHeight() == height)
				.filter(e -> e.childNode.segmentIndices != null)
				.filter(e -> e.childNode.segmentIndices[segment] == sourceNodeIndex)
				.findFirst();
	}
	
	private boolean isCoalNode(NetworkEdge e) {
		NetworkEdge sibling = getSisterEdge(e);
		for (int i = e.hasSegments.nextSetBit(0); i >= 0; i = e.hasSegments.nextSetBit(i + 1)) {
            if (sibling.hasSegments.get(i)) {
                return true;
            }
        }
		return false;
	}
	
	private NetworkEdge getSegmentNodeIndex(int segment, NetworkEdge sourceEdge) {
		if (sourceEdge.childNode.isLeaf()) {
			return sourceEdge;
		} else if (sourceEdge.childNode.isReassortment()) {
			return getSegmentNodeIndex(segment, sourceEdge.childNode.getChildEdges().get(0));
		} else {
			if (sourceEdge.childNode.getChildEdges().get(0).hasSegments.get(segment)
					&& sourceEdge.childNode.getChildEdges().get(1).hasSegments.get(segment)) {
				return sourceEdge;
			} else if (sourceEdge.childNode.getChildEdges().get(0).hasSegments.get(segment)) {
				return getSegmentNodeIndex(segment, sourceEdge.childNode.getChildEdges().get(0));
			} else {
				return getSegmentNodeIndex(segment, sourceEdge.childNode.getChildEdges().get(1));
			}
		}
	}


    private double getDelta(double edgeLength) {
        double effectiveSize = size * edgeLength;
        
        if (!gaussianInput.get()) {
            return (Randomizer.nextDouble() * effectiveSize) - (effectiveSize / 2.0);
        } else {
            return Randomizer.nextGaussian() * effectiveSize;
        }
    }


    private double getLogDeltaProb(double delta, double edgeLength) {
        double effectiveSize = size * edgeLength;
        
        if (!gaussianInput.get()) {
            if (Math.abs(delta) <= effectiveSize / 2.0) {
                return -Math.log(effectiveSize);
            } else {
                return Double.NEGATIVE_INFINITY;
            }
        } else {
            double variance = effectiveSize * effectiveSize;
            return Math.log(2) -0.5 * Math.log(2 * Math.PI * variance) - (delta * delta) / (2 * variance);
        }
    }
    
    protected BitSet getRandomConditionedSubset(BitSet sourceSegments, int requiredSegment) {
        if (sourceSegments.cardinality() < 1) {
            return null;
        }
        
        if (!sourceSegments.get(requiredSegment)) {
            return null; // Required segment not available
        }
        
        
    	BitSet destSegments = new BitSet();
        do {
        	destSegments.clear();
        
	        // Always include the required segment
	        destSegments.set(requiredSegment);
	        
	        // For all other segments, include with probability 0.5
	        for (int segIdx = sourceSegments.nextSetBit(0); segIdx != -1;
	             segIdx = sourceSegments.nextSetBit(segIdx + 1)) {
	            if (segIdx != requiredSegment && Randomizer.nextBoolean()) {
	                destSegments.set(segIdx);
	            }
	        }
        } while(destSegments.cardinality() == sourceSegments.cardinality());
        
        return destSegments;
    }
    


    protected double getLogConditionedSubsetProb(BitSet sourceSegments) {
        if (sourceSegments.cardinality() < 1) {
            return Double.NEGATIVE_INFINITY;
        }
        
        // Required segment always included (prob = 1)
        // Other segments each included with prob = 0.5
        int otherSegments = sourceSegments.cardinality() - 1;
        return otherSegments * Math.log(0.5) - Math.log(1.0 - Math.pow(0.5, otherSegments ));
    }

}
