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
public class SubTreeSlideOnNetwork extends DivertSegmentOperator {

    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<Boolean> gaussianInput = new Input<>("gaussian", "Gaussian (=true=default) or uniform delta", true);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, " +
            "i.e. -1. (when positive, gets multiplied by network-height/log2(n-taxa).", -1.0);
    

    // shadows size
    double size;
    private double limit;

	@Override
	public void initAndValidate() {
		super.initAndValidate();

        size = sizeInput.get();
        limit = limitInput.get();
	}

	@Override
	public double networkProposal() {
		
		if (divertOneSegmentInput.get()) {
			return singleSegmentProposal(); // no segments, no proposal
		}else {
			return allSegmentProposal();
		}

	}
	
	private double singleSegmentProposal() {
		double logHR = 0.0;
        // pick a random segment
        int segment = Randomizer.nextInt(network.getSegmentCount());
        
		List<NetworkEdge> potentialEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.childNode.isReassortment())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> e.hasSegments.get(segment))				
				.filter(e -> getSisterEdge(e).hasSegments.get(segment))
				.filter(e -> e.hasSegments.cardinality()==1)
			    .filter(e -> !e.parentNode.getParentEdges().get(0).isRootEdge()) // Add this
				.collect(Collectors.toList());
		
		if (potentialEdges.isEmpty()) {
			return Double.NEGATIVE_INFINITY; // no valid edges
		}
		NetworkEdge iEdge = potentialEdges.get(Randomizer.nextInt(potentialEdges.size()));
		logHR -= Math.log(1.0/potentialEdges.size());
//		System.out.println(network.getExtendedNewickVerbose(segment));
//		System.out.println(segmentTrees.get(segment)+";");
		
		// check if the parent corresponds to the segmentreeroot
//		double segmentRootHeight = getSegmentRootHeight(network.getRootEdge(), segment);
//		if (iEdge.parentNode.getHeight() == segmentRootHeight) {
//			return Double.NEGATIVE_INFINITY; // no valid edges
//		}
//		if (iEdge.parentNode.getParentEdges().get(0).isRootEdge()) {
//            return Double.NEGATIVE_INFINITY; // no valid edges
//        }
		
		NetworkEdge iEdgeChild = iEdge.childNode.getChildEdges().get(0);
		NetworkEdge jEdge = getSisterEdge(iEdge);
		
		double maxTimeSource  = Math.min(iEdge.parentNode.getHeight(), getSpouseEdge(iEdge).parentNode.getHeight());
		
	
		List<NetworkEdge> sourceEdges = collectDownEdges(iEdge, segment);
		// remove iEdge
		sourceEdges.remove(0);
		List<NetworkEdge> targetEdges = collectDownEdges(jEdge, segment);
		double minTimeDest = targetEdges.get(targetEdges.size()-1).childNode.getHeight();
		targetEdges.addAll(collectUpEdges(iEdge.parentNode.getParentEdges().get(0), segment));
		
		double minTimeSource = sourceEdges.get(sourceEdges.size()-1).childNode.getHeight();
		double maxTime = targetEdges.get(targetEdges.size()-1).childNode.getHeight();
		double minTime = Math.max(minTimeSource, minTimeDest);
	
		double destTime = Randomizer.nextDouble() * (maxTime - minTime) + minTime;

		logHR -= removeLocalReassortmentEdge(iEdge);	
		
		NetworkEdge destEdge = targetEdges.stream()				
				.filter(e -> e.childNode!=null)
				.filter(e -> e.childNode.getHeight() <= destTime)
				.filter(e -> e.parentNode.getHeight() >= destTime)
				.collect(Collectors.toList()).get(0);

		
		sourceEdges.addAll(collectUpEdges(iEdgeChild, segment, destTime));
		
		double maxStartTime = Math.min(destTime, maxTimeSource);
		
		double startTime = Randomizer.nextDouble() * (maxStartTime - minTimeSource) + minTimeSource;
		

		NetworkEdge startEdge = sourceEdges.stream()
				.filter(e -> e.childNode.getHeight() <= startTime)
				.filter(e -> e.parentNode.getHeight() >= startTime)
				.collect(Collectors.toList()).get(0);

		
		logHR += addLocalReassortmentEdge(startEdge, startTime, destEdge, destTime, segment);
//		System.out.println(segmentTrees.get(segment)+";");
		List<NetworkEdge> reversePotentialEdges = networkEdges.stream()
			    .filter(e -> !e.isRootEdge())
			    .filter(e -> e.childNode.isReassortment())
			    .filter(e -> e.parentNode.isCoalescence())
			    .filter(e -> e.hasSegments.get(segment))				
			    .filter(e -> getSisterEdge(e).hasSegments.get(segment))
			    .filter(e -> e.hasSegments.cardinality()==1)
//			    .filter(e -> e.parentNode.getHeight() != getSegmentRootHeight(network.getRootEdge(), segment)) // Add this
			    .filter(e -> !e.parentNode.getParentEdges().get(0).isRootEdge()) // Add this
			    .collect(Collectors.toList());
		
		logHR += Math.log(1.0/reversePotentialEdges.size());

		return logHR;

	}
	
	private double getSegmentRootHeight(NetworkEdge e, int segment) {
		if (e.childNode.isCoalescence()) {
			// iff both children have segment, return the height of the child
			if (e.childNode.getChildEdges().get(0).hasSegments.get(segment)
					&& e.childNode.getChildEdges().get(1).hasSegments.get(segment)) {
				return e.childNode.getHeight();
			}
			if (e.childNode.getChildEdges().get(0).hasSegments.get(segment)) {
				return getSegmentRootHeight(e.childNode.getChildEdges().get(0), segment);
			} else {
				return getSegmentRootHeight(e.childNode.getChildEdges().get(1), segment);
			}
		}else {
			return getSegmentRootHeight(e.childNode.getChildEdges().get(0), segment);
		}
	}

	private List<NetworkEdge> collectDownEdges(NetworkEdge edge, int segment) {
		List<NetworkEdge> sourceEdges = new ArrayList<>();
		
		sourceEdges.add(edge);
		
		if (edge.childNode.isLeaf()) {
			return sourceEdges;
		} else if (edge.childNode.isCoalescence()) {
			// if both children have segment, return source edges
			if (edge.childNode.getChildEdges().get(0).hasSegments.get(segment)
					&& edge.childNode.getChildEdges().get(1).hasSegments.get(segment)) {
				return sourceEdges;
			}
			// if only one child has segment, collect source edges from that child
			if (edge.childNode.getChildEdges().get(0).hasSegments.get(segment)) {
				sourceEdges.addAll(collectDownEdges(edge.childNode.getChildEdges().get(0), segment));
			} else {
				sourceEdges.addAll(collectDownEdges(edge.childNode.getChildEdges().get(1), segment));
			}
			return sourceEdges;

		} else {
			sourceEdges.addAll(collectDownEdges(edge.childNode.getChildEdges().get(0), segment));
			return sourceEdges;
		}
	}

	private List<NetworkEdge> collectUpEdges(NetworkEdge edge, int segment) {
		List<NetworkEdge> sourceEdges = new ArrayList<>();
		
		
		if (edge.isRootEdge()) {
			return sourceEdges;
		} else if (edge.parentNode.isCoalescence()) {
			sourceEdges.add(edge);
			// if both children have segment, return source edges
			if (edge.parentNode.getChildEdges().get(0).hasSegments.get(segment)
					&& edge.parentNode.getChildEdges().get(1).hasSegments.get(segment)) {
				return sourceEdges;
			}
			sourceEdges.addAll(collectUpEdges(edge.parentNode.getParentEdges().get(0), segment));
			return sourceEdges;

		} else {
			sourceEdges.add(edge);
			// if only one child has segment, collect source edges from that child
			if (edge.parentNode.getParentEdges().get(0).hasSegments.get(segment)) {
				sourceEdges.addAll(collectUpEdges(edge.parentNode.getParentEdges().get(0), segment));
			} else {
				sourceEdges.addAll(collectUpEdges(edge.parentNode.getParentEdges().get(1), segment));
			}
			return sourceEdges;
		}
	}

	private List<NetworkEdge> collectUpEdges(NetworkEdge edge, int segment, double maxHeight) {
		List<NetworkEdge> sourceEdges = new ArrayList<>();
		
		
		if (edge.isRootEdge() || edge.childNode.getHeight() > maxHeight) {
			return sourceEdges;
		} else if (edge.parentNode.isCoalescence()) {
			sourceEdges.add(edge);
			// if both children have segment, return source edges
			if (edge.parentNode.getChildEdges().get(0).hasSegments.get(segment)
					&& edge.parentNode.getChildEdges().get(1).hasSegments.get(segment)) {
				return sourceEdges;
			}
			sourceEdges.addAll(collectUpEdges(edge.parentNode.getParentEdges().get(0), segment));
			return sourceEdges;

		} else {
			sourceEdges.add(edge);
			// if only one child has segment, collect source edges from that child
			if (edge.parentNode.getParentEdges().get(0).hasSegments.get(segment)) {
				sourceEdges.addAll(collectUpEdges(edge.parentNode.getParentEdges().get(0), segment));
			} else {
				sourceEdges.addAll(collectUpEdges(edge.parentNode.getParentEdges().get(1), segment));
			}
			return sourceEdges;
		}
	}



	
	private double allSegmentProposal() {
		double logHR = 0.0;
	
        
        // pick a random segment
        int segment = Randomizer.nextInt(network.getSegmentCount());
        
		List<NetworkEdge> potentialEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> e.hasSegments.get(segment))
				.filter(e -> getSisterEdge(e).hasSegments.get(segment))
				.collect(Collectors.toList());

		NetworkEdge iEdge = potentialEdges.get(Randomizer.nextInt(potentialEdges.size()));
		
		NetworkNode iParent = iEdge.parentNode;
		NetworkNode iChild = iEdge.childNode;
		
		double sumEdgeLengths = 0.0;
		


		HashMap<NetworkEdge, Double> targetEdges = new HashMap<>();
		double minHeight = iChild.getHeight();
		double currDist = iParent.getHeight() - iChild.getHeight();
		
		final double delta = Math.exp(Math.abs(getDelta()));
		
        double deltaUp = currDist*delta;
        double deltaDown = currDist/delta;
        		
        
        
		// calculate the max height of the new child
		double maxHeight = Double.POSITIVE_INFINITY;
		
    	targetEdges.putAll(getTargetEdgesUp(iParent.getParentEdges().get(0), deltaUp, maxHeight, minHeight, segment));
    	// get the other child edge that is not iEdge
    	targetEdges.putAll(getTargetEdgesDown(getSisterEdge(iEdge), deltaDown, maxHeight, minHeight, segment));
    	
    	logHR -= Math.log(1.0/targetEdges.size());
      				        	
        List<NetworkEdge> potentialNewTargets = targetEdges.keySet().stream().collect(Collectors.toList());        
		NetworkEdge jEdge = potentialNewTargets.get(Randomizer.nextInt(potentialNewTargets.size()));
		// check if jEdge is either the same as iParent.getParentEdges().get(0) or getSisterEdge(iEdge)  		
			
		
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
		
		
		Map<NetworkEdge, Double> reverseTargetEdges = new HashMap<>();
		reverseTargetEdges.putAll(getTargetEdgesUp(iEdge.parentNode.getParentEdges().get(0), deltaUp, maxHeight, minHeight, segment));
		reverseTargetEdges.putAll(getTargetEdgesDown(getSisterEdge(iEdge), deltaDown, maxHeight, minHeight, segment));  
		
		
		logHR += Math.log(1.0/reverseTargetEdges.size());  
	    
		return logHR;
	}
	
	 
    private Map<NetworkEdge, Double> getTargetEdgesUp(NetworkEdge edge, double delta, double maxHeight, double minHeight, int segment) {
    	// initalize map
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();
        if (!edge.hasSegments.get(segment))
        	return targetEdges; // no valid targets
        
  		if (edge.isRootEdge()) {
			double newHeight = edge.childNode.getHeight() + delta;
			if (newHeight > maxHeight)
				return targetEdges; // no valid targets

			targetEdges.put(edge, newHeight);			
			return targetEdges;
		}        
		
        // get the difference between the current height and the new height
        delta -= edge.getLength();
		if (delta <= 0.0) {
			double newHeight = edge.childNode.getHeight() + delta + edge.getLength();
			if (newHeight > maxHeight)
				return targetEdges; // no valid targets
			targetEdges.put(edge, newHeight);			
			return targetEdges; 
		} else {
			// proceed to the next event
			if (edge.parentNode.isReassortment()) {
				// if reassortment, follow both parents
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(0), delta, maxHeight, minHeight, segment));
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(1), delta, maxHeight, minHeight, segment));				
			} else {
				// else it is a coalescent event
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(0), delta, maxHeight, minHeight, segment));
//				targetEdges.putAll(getTargetEdgesDown(getSisterEdge(edge), delta, maxHeight, minHeight));
			}
		}
		return targetEdges;
	}

	
    private Map<NetworkEdge, Double> getTargetEdgesDown(NetworkEdge edge, double delta, double maxHeight, double minHeight, int segment) {
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();

        if (!edge.hasSegments.get(segment))
        	return targetEdges; // no valid targets

        
        // get the difference between the current height and the new height
        delta -= edge.getLength();
        
        if (delta <= 0.0) {
        	double newHeight = edge.parentNode.getHeight() - delta - edge.getLength();
            if (newHeight < minHeight)
            	return targetEdges; // no valid targets
            
			targetEdges.put(edge, edge.parentNode.getHeight() - (delta + edge.getLength()));
			return targetEdges; 
		} else {
			// proceed to the next event
			if (edge.childNode.isLeaf()) {
				return targetEdges; // no valid targets
			}else if (edge.childNode.isReassortment()) {
				// if reassortment, follow the other parent up and the child down
//				targetEdges.putAll(getTargetEdgesUp(getSpouseEdge(edge), delta, maxHeight, minHeight));
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(0), delta, maxHeight, minHeight, segment));
			} else {
				// if coalescent event, follow both children down
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(0), delta, maxHeight, minHeight, segment));
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(1), delta, maxHeight, minHeight, segment));			
			}
		}
        return targetEdges;
	}

	double addLocalReassortmentEdge(NetworkEdge sourceEdge, double sourceTime, NetworkEdge destEdge, double destTime, int segment) {

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
		
		BitSet segsToDivert = new BitSet();

		// ensure this segment gets diverted
		segsToDivert.set(segment);
		logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
		

		networkEdges.add(reassortmentEdge);
		networkEdges.add(newEdge1);
		networkEdges.add(newEdge2);

		return logHR;
	}
	
	double removeLocalReassortmentEdge(NetworkEdge edgeToRemove) {
		double logHR = 0.0;

		NetworkNode nodeToRemove = edgeToRemove.childNode;
		NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
		NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;
		networkEdges.remove(edgeToRemove);
		networkEdges.remove(edgeToRemoveSpouse);

		// Divert segments away from chosen edge
		BitSet segsToDivert = (BitSet) edgeToRemove.hasSegments.clone();
		logHR += divertSegments(edgeToRemoveSpouse, edgeToRemove, segsToDivert);


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

//        if (!networkTerminatesAtMRCA())
//            return Double.NEGATIVE_INFINITY;

		return logHR;
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
