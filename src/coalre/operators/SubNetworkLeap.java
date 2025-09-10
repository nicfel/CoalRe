package coalre.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.util.Randomizer;
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
 * Implements the subnet slide move. General workflow:
 * 1. Choose an edge to move and a child it will carry
 * 2. Make a copy of subnet starting with this child edge
 * 3. Attach a new coppy to the new parent with randomly drawn height
 * 4. Rearrange segments
 * 5. Delete subnet starting at the child in the old position
 */
@Description("Moves the height of an internal node along the branch. " +
        "If it moves up, it can exceed the root and become a new root. " +
        "If it moves down, it may need to make a choice which branch to " +
        "slide down into.")
public class SubNetworkLeap extends DivertSegmentOperator {

    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<Boolean> gaussianInput = new Input<>("gaussian", "Gaussian (=true=default) or uniform delta", true);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, " +
            "i.e. -1. (when positive, gets multiplied by network-height/log2(n-taxa).", -1.0);
    
	public Input<Boolean> randomlySampleAttachmentEdgeInput = new Input<>("randomlySampleAttachmentEdge",
			"Randomly sample edge to attach to", true);

    // shadows size
    double size;
    private double limit;
    int maxcount;

	@Override
	public void initAndValidate() {
		super.initAndValidate();

        size = sizeInput.get();
        limit = limitInput.get();
	}
	
	@Override
	public double networkProposal() {
		try {
			if (randomlySampleAttachmentEdgeInput.get()) {
				return networkProposalRandom();
			} else {
//				System.out.println("w");
				return networkProposalWeighted();
			}
		} catch (Exception e) {
			System.err.println(e.getMessage());
			return Double.NEGATIVE_INFINITY;
		}
	}

	public double networkProposalRandom() {

		double logHR = 0.0;
		
        logHR = 0.0;
        

        // randomly sample an edge to attach to
        List<NetworkEdge> edges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> e.hasSegments.cardinality() >= 1)		
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
		
		
		
        double delta = Math.abs(getDelta(iEdge.getLength()));
//		double delta = Math.abs(getDelta(1.0));
        
        double logForwardDeltaProb = getLogDeltaProb(delta, iEdge.getLength());
//		double logForwardDeltaProb = getLogDeltaProb(delta, 1.0);

        
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
						segmentTrees.get(i).getNode(iParent.segmentIndices[i]).setHeight(targetEdges.get(jEdge));
					}	
				}
			}
		}else {
    		// keep the current segments of edge i 
    		final BitSet iSegs = (BitSet) iEdge.hasSegments.clone();
    		//topology will change
    		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
    		getTreeNodesDown(iEdge, iSegs, treeChildNodeList);
    		// remove iSegs from iEdge
    		logHR += removeSegmentsFromAncestors(iEdge, iSegs);	
    		
    		// check if parent edge is the root edge
			NetworkEdge parentEdge = iParent.getParentEdges().get(0);
    		if (parentEdge.isRootEdge()) {
	    		NetworkEdge sibEdge = getSisterEdge(iEdge);
	    		iParent.removeChildEdge(sibEdge);
	    		network.setRootEdge(sibEdge);
	    		
	    		NetworkNode newGrandParent = jEdge.parentNode;
	    		newGrandParent.removeChildEdge(jEdge);
	    		newGrandParent.addChildEdge(parentEdge);
	    		iParent.addChildEdge(jEdge);
	    		parentEdge.hasSegments = (BitSet) jEdge.hasSegments.clone();
    		}else {
    			NetworkNode grandParent = parentEdge.parentNode;	    				    			
	    		NetworkEdge sibEdge = getSisterEdge(iEdge);
	    		
	    		// reattach the sibling edge to the grandparent
	    		grandParent.removeChildEdge(parentEdge);
	    		
	    		iParent.removeChildEdge(sibEdge);		    		
	    		grandParent.addChildEdge(sibEdge);
	    		
	    		// ensure that the parent Edge has the same segments as jEdge
	    		parentEdge.hasSegments = (BitSet) jEdge.hasSegments.clone();
	    		 
	    		// readd to jEdge
	    		if (!jEdge.isRootEdge()) {
		    		NetworkNode newGrandParent = jEdge.parentNode;
		    		newGrandParent.removeChildEdge(jEdge);
		    		newGrandParent.addChildEdge(parentEdge);
		    		iParent.addChildEdge(jEdge);	
	    		}else {
		    		iParent.addChildEdge(jEdge);
		    		network.setRootEdge(parentEdge);
	    		}
    		}  
    		iEdge.parentNode.setHeight(targetEdges.get(jEdge));
    		logHR -= addSegmentsToAncestors(iEdge, iSegs);
    		
    		
    		if (reconnectSegmentTrees(treeChildNodeList, iEdge, iSegs))
    			return Double.NEGATIVE_INFINITY;
    		
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

		
			
	    double logReverseDeltaProb = getLogDeltaProb(delta, iEdge.getLength());
//		double logReverseDeltaProb = getLogDeltaProb(delta, 1.0);
	    
	    // Delta probabilities cancel out if using same edge length
	    logHR += logReverseDeltaProb - logForwardDeltaProb; // This equals 0 for symmetric case
	    	    
		int sizeAfter = (int) networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> e.hasSegments.cardinality() >= 1)	
				.filter(e -> isCoalNode(e))
				.count();

		logHR -= Math.log(1.0/sizeBefore);  // Forward
		logHR += Math.log(1.0/sizeAfter);  // Reverse

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

//        double delta = Math.abs(getDelta(iEdge.getLength()));
//        double logForwardDeltaProb = getLogDeltaProb(delta, iEdge.getLength());

        double delta = Math.abs(getDelta(1.0));
        double logForwardDeltaProb = getLogDeltaProb(delta, 1.0);

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
    		//topology will change
    		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
    		getTreeNodesDown(iEdge, iSegs, treeChildNodeList);
    		// remove iSegs from iEdge
    		logHR += removeSegmentsFromAncestors(iEdge, iSegs);	
    		
    		// check if parent edge is the root edge
			NetworkEdge parentEdge = iParent.getParentEdges().get(0);
    		if (parentEdge.isRootEdge()) {
	    		NetworkEdge sibEdge = getSisterEdge(iEdge);
	    		iParent.removeChildEdge(sibEdge);
	    		network.setRootEdge(sibEdge);
	    		
	    		NetworkNode newGrandParent = jEdge.parentNode;
	    		newGrandParent.removeChildEdge(jEdge);
	    		newGrandParent.addChildEdge(parentEdge);
	    		iParent.addChildEdge(jEdge);
	    		parentEdge.hasSegments = (BitSet) jEdge.hasSegments.clone();
    		}else {
    			NetworkNode grandParent = parentEdge.parentNode;	    				    			
	    		NetworkEdge sibEdge = getSisterEdge(iEdge);
	    		
	    		// reattach the sibling edge to the grandparent
	    		grandParent.removeChildEdge(parentEdge);
	    		
	    		iParent.removeChildEdge(sibEdge);		    		
	    		grandParent.addChildEdge(sibEdge);
	    		
	    		// ensure that the parent Edge has the same segments as jEdge
	    		parentEdge.hasSegments = (BitSet) jEdge.hasSegments.clone();
	    		 
	    		// readd to jEdge
	    		if (!jEdge.isRootEdge()) {
		    		NetworkNode newGrandParent = jEdge.parentNode;
		    		newGrandParent.removeChildEdge(jEdge);
		    		newGrandParent.addChildEdge(parentEdge);
		    		iParent.addChildEdge(jEdge);	
	    		}else {
		    		iParent.addChildEdge(jEdge);
		    		network.setRootEdge(parentEdge);
	    		}
    		}  
    		iEdge.parentNode.setHeight(targetEdges.get(jEdge));
    		logHR -= addSegmentsToAncestors(iEdge, iSegs);
    		
    		
    		if (reconnectSegmentTrees(treeChildNodeList, iEdge, iSegs))
    			return Double.NEGATIVE_INFINITY;
    		
		}
		maxcount = 0; // should not be necessary TODO: investigate why this is needed

		Map<NetworkEdge, Double> reverseTargetEdges = new HashMap<>();
		List<NetworkEdge> visitedNodes2 = new ArrayList<>();
		reverseTargetEdges.putAll(getTargetEdgesUp(iEdge.parentNode.getParentEdges().get(0), delta, minHeight, iEdge, segment));
		reverseTargetEdges.putAll(getTargetEdgesDown(getSisterEdge(iEdge), delta, minHeight, iEdge, segment));  
		
				
		logHR += Math.log(1.0/reverseTargetEdges.size());  
		
			
//	    double logReverseDeltaProb = getLogDeltaProb(delta, iEdge.getLength());
	    double logReverseDeltaProb = getLogDeltaProb(delta, 1.0);
	    
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
//        System.out.println("a " + edge.childNode.getHeight() + "  " + edge .parentNode.getHeight() + " " + delta + " " + ii);
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
