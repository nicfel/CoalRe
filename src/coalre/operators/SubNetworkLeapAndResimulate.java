package coalre.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.util.Randomizer;
import coalre.distribution.NetworkEvent;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;


/**
 * Implements the subnet slide move with resimulation. General workflow:
 * 1. Choose an edge to move and a child it will carry
 * 2. Make a copy of subnet starting with this child edge
 * 3. Attach a new copy to the new parent with randomly drawn height
 * 4. Resimulate segments forward using coalescent with reassortment
 * 5. Delete subnet starting at the child in the old position
 */
@Description("Moves the height of an internal node along the branch using resimulation. " +
        "If it moves up, it can exceed the root and become a new root. " +
        "If it moves down, it may need to make a choice which branch to " +
        "slide down into.")
public class SubNetworkLeapAndResimulate extends DivertSegmentAndResimulate {

    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<Boolean> gaussianInput = new Input<>("gaussian", "Gaussian (=true=default) or uniform delta", true);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, " +
            "i.e. -1. (when positive, gets multiplied by network-height/log2(n-taxa).", -1.0);
    
	public Input<Boolean> randomlySampleAttachmentEdgeInput = new Input<>("randomlySampleAttachmentEdge",
			"Randomly sample edge to attach to", true);
	public Input<Boolean> useEdgeLengthInput = new Input<>("useEdgeLength",
			"Use edge length to scale the size of the move", false);

    // shadows size
    double size;
    private double limit;
    int maxcount;
    
    boolean useSingleSegmentDiversion = false;

	@Override
	public void initAndValidate() {
		super.initAndValidate();

        size = sizeInput.get();
        limit = limitInput.get();
        
	}
	
	@Override
	public double networkProposal() {
		double logHR = 0.0;
//		try {
			if (randomlySampleAttachmentEdgeInput.get()) {
				logHR = networkProposalRandom();
			} else {
				logHR = networkProposalWeighted();
			}
//		} catch (Exception e) {
//			System.err.println(e.getMessage());
//			return Double.NEGATIVE_INFINITY;
//		}
			

		return logHR;

	}

	public double networkProposalRandom() {

		double logHR = 0.0;
		
        logHR = 0.0;
        

        // randomly sample an edge to attach to
        List<NetworkEdge> edges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> isCoalNode(e))
				.collect(Collectors.toList());

		NetworkEdge iEdge = edges.get(Randomizer.nextInt(edges.size()));

		int sizeBefore = edges.size();
		
		NetworkNode iParent = iEdge.parentNode;		
		NetworkNode iChild = iEdge.childNode;
		
		BitSet coalSegments = (BitSet) iEdge.hasSegments.clone();
		coalSegments.and(getSisterEdge(iEdge).hasSegments);
//		System.out.println(coalSegments);
		// pick a random segment on the source edge
		int segment = Randomizer.nextInt(network.getSegmentCount());
		while (!coalSegments.get(segment)) {
			segment = Randomizer.nextInt(network.getSegmentCount());
		}
		logHR -= Math.log(1.0 / coalSegments.cardinality()); // pick random segment
		
		double delta;
		double logForwardDeltaProb;
		if (useEdgeLengthInput.get()) {
	        delta = Math.abs(getDelta(iEdge.getLength()));        
	        logForwardDeltaProb = getLogDeltaProb(delta, iEdge.getLength());
		}else {
			delta = Math.abs(getDelta(1.0));
			logForwardDeltaProb = getLogDeltaProb(delta, 1.0);
		}

        
        // get all potential reattachment Edges
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();

		// calculate the max height of the new child
		double minHeight = iChild.getHeight();

		maxcount = 0; // should not be neccessary TODO: investigate why this is needed
		
//		System.out.println("a");
    	targetEdges.putAll(getTargetEdgesUp(iParent.getParentEdges().get(0), delta, minHeight, iEdge, segment));
//    	System.out.println("b");
    	// get the other child edge that is not iEdge
    	targetEdges.putAll(getTargetEdgesDown(getSisterEdge(iEdge), delta, minHeight, iEdge, segment));
    	logHR -= Math.log(1.0/targetEdges.size());
      				        	
        List<NetworkEdge> potentialNewTargets = targetEdges.keySet().stream().collect(Collectors.toList());        
		NetworkEdge jEdge = potentialNewTargets.get(Randomizer.nextInt(potentialNewTargets.size()));
		if (jEdge == iParent.getParentEdges().get(0) || jEdge == getSisterEdge(iEdge)) {
			// change the height of jParent, but no topological changes
			iParent.setHeight(targetEdges.get(jEdge));
			if (iParent.segmentIndices != null && iParent.isCoalescence()) {
				for (int i = 0; i < iParent.segmentIndices.length; i++) {
					if (iParent.getChildEdges().get(0).hasSegments.get(i)
							&& iParent.getChildEdges().get(1).hasSegments.get(i)) {
//						System.out.println(iParent.getHeight());
						segmentTrees.get(i).getNode(iParent.segmentIndices[i]).setHeight(targetEdges.get(jEdge));
					}	
				}
			}
		}else {
    		// keep the current segments of edge i 
    		final BitSet iSegs = (BitSet) iEdge.hasSegments.clone();
    		// Use resimulation approach instead of simple add/remove
    		logHR += divertSegmentsWithResimulation(iEdge, jEdge, iSegs, targetEdges.get(jEdge));
    		
		}
		
		maxcount = 0; // should not be necessary TODO: investigate why this is needed
		Map<NetworkEdge, Double> reverseTargetEdges = new HashMap<>();
		List<NetworkEdge> visitedNodes2 = new ArrayList<>();
//		System.out.println("c");
		reverseTargetEdges.putAll(getTargetEdgesUp(iEdge.parentNode.getParentEdges().get(0), delta, minHeight, iEdge, segment));
//		System.out.println("d");
		reverseTargetEdges.putAll(getTargetEdgesDown(getSisterEdge(iEdge), delta, minHeight, iEdge, segment));  
		
				
		logHR += Math.log(1.0/reverseTargetEdges.size());  
		// compute the reverse number of segments
		BitSet coalSegmentsReverse = (BitSet) iEdge.hasSegments.clone();
		coalSegmentsReverse.and(getSisterEdge(iEdge).hasSegments);
		
		logHR += Math.log(1.0 / coalSegmentsReverse.cardinality()); // pick random segment

		
		double logReverseDeltaProb;
		if (useEdgeLengthInput.get())
			logReverseDeltaProb = getLogDeltaProb(delta, iEdge.getLength());
		else
	    	logReverseDeltaProb = getLogDeltaProb(delta, 1.0);

	    // Delta probabilities cancel out if using same edge length
	    logHR += logReverseDeltaProb - logForwardDeltaProb; // This equals 0 for symmetric case
	    	    
		int sizeAfter = (int) networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> isCoalNode(e))
				.count();

		logHR -= Math.log(1.0/sizeBefore);  // Forward
		logHR += Math.log(1.0/sizeAfter);  // Reverse
							 
//		if (sizeAfter != sizeBefore)
//			System.out.println(sizeBefore - sizeAfter);
//		System.out.println("");
		return logHR;
	}
	
	public double networkProposalWeighted() {

		double logHR = 0.0;
		
        logHR = 0.0;
        
        List<Integer> validEdges = new ArrayList<>();
        List<Double> edgeLength = new ArrayList<>();
        double sumLength = 0.0;

		for (int i = 0; i < networkEdges.size(); i++ ) {
			if (!networkEdges.get(i).isRootEdge() 
					&& networkEdges.get(i).parentNode.isCoalescence()
					&& networkEdges.get(i).hasSegments.cardinality() >= 1
					&& isCoalNode(networkEdges.get(i))) {
				validEdges.add(i);
				edgeLength.add(networkEdges.get(i).getLength());
				sumLength += networkEdges.get(i).getLength();
			}
		}
		
		// sample proportionally to the edge length
		double r = Randomizer.nextDouble() * sumLength;
		int idx = 0;
		for (int i = 0; i < edgeLength.size(); i++) {
			r -= edgeLength.get(i);
			if (r <= 0) {
				idx = i;
				break;
			}
		
		}
		        
		logHR -= Math.log(edgeLength.get(idx)/sumLength);  // Forward
        

		NetworkEdge iEdge = networkEdges.get(validEdges.get(idx));

		
		NetworkNode iParent = iEdge.parentNode;		
		NetworkNode iChild = iEdge.childNode;
		
		
		BitSet coalSegments = (BitSet) iEdge.hasSegments.clone();
		coalSegments.and(getSisterEdge(iEdge).hasSegments);
		// pick a random segment on the source edge
		int segment = Randomizer.nextInt(network.getSegmentCount());
		while (!coalSegments.get(segment)) {
			segment = Randomizer.nextInt(network.getSegmentCount());
		}
		logHR -= Math.log(1.0 / coalSegments.cardinality()); // pick random segment

		double delta;
		double logForwardDeltaProb;
		if (useEdgeLengthInput.get()) {
	        delta = Math.abs(getDelta(iEdge.getLength()));        
	        logForwardDeltaProb = getLogDeltaProb(delta, iEdge.getLength());
		}else {
			delta = Math.abs(getDelta(1.0));
			logForwardDeltaProb = getLogDeltaProb(delta, 1.0);
		}

        
        // get all potential reattachment Edges
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();

		// calculate the max height of the new child
		double minHeight = iChild.getHeight();
		maxcount = 0; // should not be necessary TODO: investigate why this is needed

		List<NetworkEdge> visitedNodes = new ArrayList<>();
    	targetEdges.putAll(getTargetEdgesUp(iParent.getParentEdges().get(0), delta, minHeight, iEdge, segment));
    	// get the other child edge that is not iEdge
    	targetEdges.putAll(getTargetEdgesDown(getSisterEdge(iEdge), delta, minHeight, iEdge, segment));
    	logHR -= Math.log(1.0/targetEdges.size());
      				        	
        List<NetworkEdge> potentialNewTargets = targetEdges.keySet().stream().collect(Collectors.toList());        
		NetworkEdge jEdge = potentialNewTargets.get(Randomizer.nextInt(potentialNewTargets.size()));
		
		if (jEdge == iParent.getParentEdges().get(0) || jEdge == getSisterEdge(iEdge)) {
			// change the height of jParent, but no topological changes
			iParent.setHeight(targetEdges.get(jEdge));
			if (iParent.segmentIndices != null && iParent.isCoalescence()) {
				for (int i = 0; i < iParent.segmentIndices.length; i++) {
					if (iParent.getChildEdges().get(0).hasSegments.get(i)
							&& iParent.getChildEdges().get(1).hasSegments.get(i)) {
						segmentTrees.get(i).getNode(iParent.segmentIndices[i]).setHeight(targetEdges.get(jEdge));
					}	
				}
			}
		}else {
    		// keep the current segments of edge i 
    		final BitSet iSegs = (BitSet) iEdge.hasSegments.clone();
    		//topology will change - use resimulation approach
    		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
    		getTreeNodesDown(iEdge, iSegs, treeChildNodeList);
    		
    		// Use resimulation approach instead of simple add/remove
    		logHR += divertSegmentsWithResimulation(iEdge, jEdge, iSegs, targetEdges.get(jEdge));
    		
		}
		maxcount = 0; // should not be necessary TODO: investigate why this is needed

		Map<NetworkEdge, Double> reverseTargetEdges = new HashMap<>();
		List<NetworkEdge> visitedNodes2 = new ArrayList<>();
		reverseTargetEdges.putAll(getTargetEdgesUp(iEdge.parentNode.getParentEdges().get(0), delta, minHeight, iEdge, segment));
		reverseTargetEdges.putAll(getTargetEdgesDown(getSisterEdge(iEdge), delta, minHeight, iEdge, segment));  
		
				
		logHR += Math.log(1.0/reverseTargetEdges.size());  
		
			
		double logReverseDeltaProb;
		if (useEdgeLengthInput.get())
			logReverseDeltaProb = getLogDeltaProb(delta, iEdge.getLength());
		else
	    	logReverseDeltaProb = getLogDeltaProb(delta, 1.0);
	    
	    // Delta probabilities cancel out if using same edge length
	    logHR += logReverseDeltaProb - logForwardDeltaProb; // This equals 0 for symmetric case
	    	    
	    
	    
		BitSet coalSegmentsReverse = (BitSet) iEdge.hasSegments.clone();
		coalSegmentsReverse.and(getSisterEdge(iEdge).hasSegments);
		
		logHR += Math.log(1.0 / coalSegmentsReverse.cardinality()); // pick random segment

	    

        double sumLengthReverse = 0.0;

		for (int i = 0; i < networkEdges.size(); i++ ) {
			if (!networkEdges.get(i).isRootEdge() 
					&& networkEdges.get(i).parentNode.isCoalescence()
					&& networkEdges.get(i).hasSegments.cardinality() >= 1
					&& isCoalNode(networkEdges.get(i))) {
				sumLengthReverse += networkEdges.get(i).getLength();
			}
		}
		
		logHR += Math.log(iEdge.getLength()/sumLengthReverse);  // Reverse
							    
		return logHR;
	}
	

	private Map<NetworkEdge, Double> getTargetEdgesUp(NetworkEdge edge, double delta, double minHeight, NetworkEdge iEdge, int segment) {

    	maxcount++;
		if (maxcount > 1000000) {
//			System.out.println(network);
//			System.out.println(iEdge.childNode.getHeight());
			throw new IllegalStateException("Too many iterations in getTargetEdgesUp: " + maxcount);
		}
    	// initalize map
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();
		if (!edge.hasSegments.get(segment)) {
			return targetEdges; // do not move the edge itself
		}

        if (edge==iEdge)
        	return targetEdges; // do not move the edge itself
//        ii++;
  		if (edge.isRootEdge()) {
			double newHeight = edge.childNode.getHeight() + delta;
			targetEdges.put(edge, newHeight);			
			return targetEdges;
		}        
  			
        // get the difference between the current height and the new height
        delta -= edge.getLength();
        if (edge.getLength()<=0) {
        	throw new IllegalArgumentException("Edge length cannot be negative: " + edge.getLength());
        }
		if (delta <= 0.0) {
			double newHeight = edge.childNode.getHeight() + delta + edge.getLength();
			targetEdges.put(edge, newHeight);			
			return targetEdges; 
		} else {
			// proceed to the next event
			if (edge.parentNode.isReassortment()) {
				// if reassortment, follow both parents
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(0), delta, minHeight, iEdge, segment));
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(1), delta, minHeight, iEdge, segment));				
			} else {
				// else it is a coalescent event
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(0), delta, minHeight, iEdge, segment));
				targetEdges.putAll(getTargetEdgesDown(getSisterEdge(edge), delta, minHeight, iEdge, segment));
			}
		}
		return targetEdges;
	}
	
    private Map<NetworkEdge, Double> getTargetEdgesDown(NetworkEdge edge, double delta, double minHeight, NetworkEdge iEdge, int segment) {
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();
//        ii++;
        
		if (!edge.hasSegments.get(segment)) {
			return targetEdges; // do not move the edge itself
		}
        
    	maxcount++;
		if (maxcount > 1000000) {
//			System.out.println(network);
//			System.out.println(iEdge.childNode.getHeight());

			throw new IllegalStateException("Too many iterations in getTargetEdgesUp: " + maxcount);
		}

		if (edge == iEdge)
			return targetEdges; // do not move the edge itself
        
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
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(0), delta, minHeight, iEdge, segment));
			} else {
				// if coalescent event, follow both children down
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(0), delta, minHeight, iEdge, segment));
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(1), delta, minHeight, iEdge, segment));			
			}
		}
        return targetEdges;
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
            return -0.5 * Math.log(2 * Math.PI * variance) - (delta * delta) / (2 * variance);
        }
    }


    private double getDelta() {
        if (!gaussianInput.get()) {
            return (Randomizer.nextDouble() * size) - (size / 2.0);
        } else {
            return Randomizer.nextGaussian() * size;
        }
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

	/**
	 * Divert segments from sourceEdge to destEdge using resimulation approach.
	 * This method handles the topological changes and resimulates the history forward.
	 */
	protected double divertSegmentsWithResimulation(NetworkEdge sourceEdge, NetworkEdge destEdge, BitSet segsToDivert, double newHeight) {
		double logHR = 0.0;
		
		coalescentDistr.intervals.update();
		// Get network events to determine lineages at different times
		List<NetworkEvent> networkEventList;

		
        // make a temporary new node on the source edge that is below the new height
        NetworkNode sourceNode = new NetworkNode();
        NetworkNode destNode = new NetworkNode();
        NetworkNode oldSourceParent = sourceEdge.parentNode;
        
        double sourceHeight = Math.min(newHeight, sourceEdge.parentNode.getHeight());
        
        sourceNode.setHeight((sourceHeight + sourceEdge.childNode.getHeight()) / 2.0);
        destNode.setHeight(newHeight);
        
        
        NetworkEdge newEdgeLeft = new NetworkEdge();
        NetworkEdge newEdgeRight = new NetworkEdge();
        NetworkEdge newDestChildEdge = new NetworkEdge();
        newDestChildEdge.hasSegments = (BitSet) destEdge.hasSegments.clone();
        newEdgeLeft.hasSegments = (BitSet) sourceEdge.hasSegments.clone();
        newEdgeRight.hasSegments = new BitSet();
        
        NetworkNode destChild = destEdge.childNode;
        destChild.removeParentEdge(destEdge);
        destChild.addParentEdge(newDestChildEdge);
        
        networkEdges.add(newEdgeLeft);
        networkEdges.add(newEdgeRight);
        networkEdges.add(newDestChildEdge);
        
        
        oldSourceParent.removeChildEdge(sourceEdge);
        oldSourceParent.addChildEdge(newEdgeLeft);

        sourceNode.addParentEdge(newEdgeLeft);
        sourceNode.addParentEdge(newEdgeRight);
        destNode.addChildEdge(newDestChildEdge);        
        destNode.addChildEdge(newEdgeRight);
        destNode.addParentEdge(destEdge);
        sourceNode.addChildEdge(sourceEdge);
        
        
        
        
		
		// Track tree nodes before removal
		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
		getTreeNodesDown(sourceEdge, segsToDivert, treeChildNodeList);

		BitSet segsToDivertCopy = (BitSet) segsToDivert.clone();

		
		if (useSingleSegmentDiversion) {
			// make a defensive implementation by selection one segment at a time to divert, choose the order randomly
			for (int i = 0; i < segsToDivert.cardinality(); i++) {
				// Get network events to determine lineages at different times
				networkEventList = coalescentDistr.intervals.getNetworkEventList();
	
						
	
				// select a random segment from segsToDivertCopy
				int segIdx = -1;
				int randIdx = Randomizer.nextInt(segsToDivertCopy.cardinality());
				for (int j = 0; j <= randIdx; j++) {
					segIdx = segsToDivertCopy.nextSetBit(segIdx + 1);
				}
				segsToDivertCopy.set(segIdx, false);
				
				BitSet currentSegsToDivert = new BitSet();
				currentSegsToDivert.set(segIdx);
				
		        // Prepare for resimulation
				BitSet segsToAdd = (BitSet) currentSegsToDivert.clone();
				List<NetworkEdge> activeEdges = new ArrayList<>();
				List<BitSet> segsToAddList = new ArrayList<>();
				
				List<NetworkEdge> inactiveEdges = new ArrayList<>();
				List<BitSet> segsToRemoveList = new ArrayList<>();
	
			
				// Set up for resimulation
				activeEdges.add(newEdgeRight);
				segsToAddList.add(segsToAdd);
				
				inactiveEdges.add(newEdgeLeft);
				segsToRemoveList.add((BitSet) currentSegsToDivert.clone());
				
				// Track the current time (start at the bottom of the lowest edge)
				double currentTime = sourceNode.getHeight();
		//		System.out.println(network);
		//        System.out.println(sourceHeight +", "+ newHeight + ", "+ sourceNode.getHeight());
		//        System.exit(0);
		
				
				// Use resimulation approach
				logHR += divertSegmentsToAncestors(activeEdges, inactiveEdges, segsToAddList, segsToRemoveList, 
						currentTime, networkEventList, true, true);
		
				if (reconnectSegmentTrees(treeChildNodeList, newEdgeRight, currentSegsToDivert))
					return Double.NEGATIVE_INFINITY;
					
	
				cleanEmptyEdgesTopDown();
				if (segsToDivert.cardinality()!=i+1) {
					coalescentDistr.intervals.update();
				}
			}
		}else {
			// Get network events to determine lineages at different times
			networkEventList = coalescentDistr.intervals.getNetworkEventList();
	        // Prepare for resimulation
			BitSet segsToAdd = (BitSet) segsToDivertCopy.clone();
			List<NetworkEdge> activeEdges = new ArrayList<>();
			List<BitSet> segsToAddList = new ArrayList<>();
			
			List<NetworkEdge> inactiveEdges = new ArrayList<>();
			List<BitSet> segsToRemoveList = new ArrayList<>();

		
			// Set up for resimulation
			activeEdges.add(newEdgeRight);
			segsToAddList.add(segsToAdd);
			
			inactiveEdges.add(newEdgeLeft);
			segsToRemoveList.add((BitSet) segsToDivertCopy.clone());
			
			// Track the current time (start at the bottom of the lowest edge)
			double currentTime = sourceNode.getHeight();
	//		System.out.println(network);
	//        System.out.println(sourceHeight +", "+ newHeight + ", "+ sourceNode.getHeight());
	//        System.exit(0);
	
			
			// Use resimulation approach
			logHR += divertSegmentsToAncestors(activeEdges, inactiveEdges, segsToAddList, segsToRemoveList, 
					currentTime, networkEventList, true, true);
	
			if (reconnectSegmentTrees(treeChildNodeList, newEdgeRight, segsToDivertCopy))
				return Double.NEGATIVE_INFINITY;
				

			logHR += cleanEmptyEdgesTopDown();
		}
		

		return logHR;
	}

    
    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(size);
            final double f = Math.exp(delta);
            if( limit > 0 ) {
                final Network network = networkInput.get();
                final double h = network.getRootEdge().childNode.getHeight();
                final double k = Math.log(network.getLeafNodes().size()) / Math.log(2);
                final double lim = (h / k) * limit;
                if( f <= lim ) {
                    size = f;
                }
            } else {
               size = f;
            }
        }

    }

    @Override
    public double getCoercableParameterValue() {
        return size;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        size = value;
    }
    
    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;

        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        final double newDelta = size * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try decreasing size to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try increasing size to about " + formatter.format(newDelta);
        } else return "";
    }
}

