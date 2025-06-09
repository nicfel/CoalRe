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
public class SubNetworkSlide extends DivertSegmentOperator {

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

		double logHR = 0.0;
		network.startEditing(this);
//		System.out.println();
		
//		System.out.println(network);
		
		// 1. Choose a random edge, avoiding root
		List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());
		
		final List<NetworkEdge> possibleEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.collect(Collectors.toList());
		
		final int possibleNodes = possibleEdges.size();
        if (possibleNodes < 1) {
            return Double.NEGATIVE_INFINITY;
        }

        logHR = 0.0;

		NetworkEdge iEdge = possibleEdges.get(Randomizer.nextInt(possibleNodes));
		NetworkNode iParent = iEdge.parentNode;
		NetworkNode iChild = iEdge.childNode;

		
        final double delta = Math.abs(getDelta());

        // get all potential reattachment Edges
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();
        if (!iEdge.isRootEdge() && iParent.isReassortment()) {

        	NetworkEdge iParentEdge = iParent.getParentEdges().get(Randomizer.nextInt(2));    
        	
    		double maxHeight = iParentEdge.parentNode.getHeight();
    		double minHeight = 0;
        	
        	// get the other pare
        	targetEdges.putAll(getTargetEdgesUp(getSpouseEdge(iParentEdge), delta, maxHeight, minHeight));
        	targetEdges.putAll(getTargetEdgesDown(iEdge, delta, maxHeight, minHeight));  
        	        	
        	// if list is empty return negative infinity
			if (targetEdges.isEmpty()) {
				return Double.NEGATIVE_INFINITY; // no valid targets
			}
        	
			logHR -= Math.log(1.0/targetEdges.size());  
        	
        	// randomly select an edge
            List<NetworkEdge> potentialNewTargets = targetEdges.keySet().stream().collect(Collectors.toList());        
    		NetworkEdge jEdge = potentialNewTargets.get(Randomizer.nextInt(potentialNewTargets.size()));
    		// keep the current segments of edge i 

			if (jEdge == iParentEdge.childNode.getChildEdges().get(0) || jEdge == getSpouseEdge(iParentEdge)) {
				// change the height of jParent
				iParent.setHeight(targetEdges.get(jEdge));
			}else {			
				
	    		// keep the current segments of edge i 
	    		final BitSet iParentSegs = (BitSet) iParentEdge.hasSegments.clone();
	    		
	            // check that jEdge has all segments in iSegs
				for (int i = 0; i < iParentSegs.length(); i++) {
					if (iParentSegs.get(i) && !jEdge.hasSegments.get(i)) {
						return Double.NEGATIVE_INFINITY; // cannot move segments that are not present in jEdge
					}
				}	

				// remove iParentEdge from child node
				NetworkEdge grandChildEdge = iParentEdge.childNode.getChildEdges().get(0);
				NetworkEdge otherParentEdge = getSpouseEdge(iParentEdge);
				NetworkNode uncle = otherParentEdge.parentNode;
				
	    		// remove iSegs from iEdge				
	    		logHR += addSegmentsToAncestors(otherParentEdge, iParentSegs);	   
				
				uncle.removeChildEdge(otherParentEdge);
				iParentEdge.childNode.removeChildEdge(grandChildEdge);
				uncle.addChildEdge(grandChildEdge);

				// read to jEdge
				NetworkNode newParent = jEdge.parentNode;
				
				newParent.removeChildEdge(jEdge);
				newParent.addChildEdge(otherParentEdge);
				iParentEdge.childNode.addChildEdge(jEdge);
				iParentEdge.childNode.setHeight(targetEdges.get(jEdge));

				otherParentEdge.hasSegments = (BitSet) jEdge.hasSegments.clone();
				logHR -= removeSegmentsFromAncestors(otherParentEdge, iParentSegs);	    
			
			}

			Map<NetworkEdge, Double> reverseTargetEdges = new HashMap<>();
            reverseTargetEdges.putAll(getTargetEdgesUp(getSpouseEdge(iParentEdge), delta, maxHeight, minHeight));
            reverseTargetEdges.putAll(getTargetEdgesDown(jEdge, delta, maxHeight, minHeight));
                        
			logHR += Math.log(1.0/reverseTargetEdges.size());  

			
			return Double.NEGATIVE_INFINITY; // cannot move segments that are not present in jEdge

        }else{
    		// calculate the max height of the new child
    		double maxHeight = Double.POSITIVE_INFINITY;
    		double minHeight = iChild.getHeight();
        	targetEdges.putAll(getTargetEdgesUp(iParent.getParentEdges().get(0), delta, maxHeight, minHeight));
        	// get the other child edge that is not iEdge
        	targetEdges.putAll(getTargetEdgesDown(getSisterEdge(iEdge), delta, maxHeight, minHeight));
        	
        	logHR -= Math.log(1.0/targetEdges.size());
          				        	
            List<NetworkEdge> potentialNewTargets = targetEdges.keySet().stream().collect(Collectors.toList());        
    		NetworkEdge jEdge = potentialNewTargets.get(Randomizer.nextInt(potentialNewTargets.size()));
    		// check if jEdge is either the same as iParent.getParentEdges().get(0) or getSisterEdge(iEdge)  		
    			
    		
			if (jEdge == iParent.getParentEdges().get(0) || jEdge == getSisterEdge(iEdge)) {
				// change the height of jParent
				iParent.setHeight(targetEdges.get(jEdge));
			}else {
	    		// keep the current segments of edge i 
	    		final BitSet iSegs = (BitSet) iEdge.hasSegments.clone();
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
			}
			
			
			Map<NetworkEdge, Double> reverseTargetEdges = new HashMap<>();
			reverseTargetEdges.putAll(getTargetEdgesUp(iEdge.parentNode.getParentEdges().get(0), delta, maxHeight, minHeight));
			reverseTargetEdges.putAll(getTargetEdgesDown(getSisterEdge(iEdge), delta, maxHeight, minHeight));  
			
			logHR += Math.log(1.0/reverseTargetEdges.size());  
			
        }
        		
	    if(!networkTerminatesAtMRCA())
	    	return Double.NEGATIVE_INFINITY;
	    
//	    // for each coalescent event, check that the parent edge is segments that is the union of the child edges
//	    for (NetworkNode n : network.getInternalNodes()) {
//	    	if (n.isCoalescence()) {
//	    		BitSet child1 = (BitSet) n.getChildEdges().get(0).hasSegments.clone();
//	    		child1.or(n.getChildEdges().get(1).hasSegments);
//	    		// check that they are the same
//				if (!child1.equals(n.getParentEdges().get(0).hasSegments)) {
//					System.out.println("Error: parent edge segments do not match child edges: "
//							+ n.getParentEdges().get(0).hasSegments + " " + child1);
//					System.exit(0);
//				}	    		
//	    	}else if (n.isReassortment()) {
//	    		BitSet parent1 = (BitSet) n.getParentEdges().get(0).hasSegments.clone();
//	    		parent1.or(n.getParentEdges().get(1).hasSegments);
//	    		// check that they are the same
//				if (!parent1.equals(n.getChildEdges().get(0).hasSegments)) {
//					System.out.println("Error: child edge segments do not match parent edges: "
//							+ n.getChildEdges().get(0).hasSegments + " " + parent1);
//					System.exit(0);
//				}
//	    	}
//	    }

		return logHR;
	}
	
	 
    private Map<NetworkEdge, Double> getTargetEdgesUp(NetworkEdge edge, double delta, double maxHeight, double minHeight) {
    	// initalize map
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();
        
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
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(0), delta, maxHeight, minHeight));
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(1), delta, maxHeight, minHeight));				
			} else {
				// else it is a coalescent event
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(0), delta, maxHeight, minHeight));
//				targetEdges.putAll(getTargetEdgesDown(getSisterEdge(edge), delta, maxHeight, minHeight));
			}
		}
		return targetEdges;
	}

	
    private Map<NetworkEdge, Double> getTargetEdgesDown(NetworkEdge edge, double delta, double maxHeight, double minHeight) {
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();
        
        
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
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(0), delta, maxHeight, minHeight));
			} else {
				// if coalescent event, follow both children down
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(0), delta, maxHeight, minHeight));
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(1), delta, maxHeight, minHeight));			
			}
		}
        return targetEdges;
	}

	public double checkAncestralSegments(NetworkEdge edge) {
    	double logHR = 0.0;
    	
    	NetworkNode parentNode = edge.parentNode;
    	if (parentNode == null) {
    		return logHR;
    	}
    	
    	if (parentNode.getParentCount() > 1) {
    		NetworkEdge parentEdge = parentNode.getParentEdges().get(0);
    		NetworkEdge spouseEdge = parentNode.getParentEdges().get(1);
    		
    		BitSet segsToAdd = (BitSet)edge.hasSegments.clone();
    		segsToAdd.andNot(parentEdge.hasSegments);
    		segsToAdd.andNot(spouseEdge.hasSegments);
    		
    		if (!segsToAdd.isEmpty()) {
    			if (Randomizer.nextBoolean()) {
    				logHR -= addSegmentsToAncestors(parentEdge, segsToAdd);
    			} else {
    				logHR -= addSegmentsToAncestors(spouseEdge, segsToAdd);
    			}
    		}
    		
    		BitSet segsToRemoveParent = (BitSet)parentEdge.hasSegments.clone();
    		segsToRemoveParent.andNot(edge.hasSegments);
    		if (!segsToRemoveParent.isEmpty()) {
    			logHR += removeSegmentsFromAncestors(parentEdge, segsToRemoveParent);
    		}
    		
    		BitSet segsToRemoveSpouse= (BitSet)spouseEdge.hasSegments.clone();
    		segsToRemoveSpouse.andNot(edge.hasSegments);
    		if (!segsToRemoveSpouse.isEmpty()) {
    			logHR += removeSegmentsFromAncestors(spouseEdge, segsToRemoveSpouse);
    		}
	
    	} else {
    		NetworkEdge parentEdge = parentNode.getParentEdges().get(0);
    		BitSet segsToAdd = (BitSet)edge.hasSegments.clone();
    		segsToAdd.andNot(parentEdge.hasSegments);
    		
    		if (!segsToAdd.isEmpty()) {
    			logHR -= addSegmentsToAncestors(parentEdge, segsToAdd);
    		}
    		
    		BitSet segsToRemoveParent = (BitSet)parentEdge.hasSegments.clone();
    		segsToRemoveParent.andNot(edge.hasSegments);
    		if (!segsToRemoveParent.isEmpty()) {
    			logHR += removeSegmentsFromAncestors(parentEdge, segsToRemoveParent);
    		}
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
