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
public class SubNetworkLeapAddRemove extends AddRemoveReassortment {

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

	@Override
	public void initAndValidate() {
		super.initAndValidate();

        size = sizeInput.get();
        limit = limitInput.get();
	}
	int ii=0;
	
	@Override
	public double networkProposal() {
		return 0.0;
//		double logHR = 0.0;
//		if (Randomizer.nextBoolean()) {
//			logHR = addReassortmentEvent();
//		} else {
//			try {
//				logHR = removeRecombination();
//			} catch (Exception e) {
//				return Double.NEGATIVE_INFINITY;
//			}
//		}		

	}

	
	public double addReassortmentEvent() {

		double logHR = 0.0;
		        
        // pick a time between root and time 0
		double sourceTime = Randomizer.nextDouble() * network.getRootEdge().childNode.getHeight();
		logHR -= Math.log(1.0/network.getRootEdge().childNode.getHeight());        
        
		List<NetworkEdge> edges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> e.childNode.getHeight() < sourceTime)
				.filter(e -> e.parentNode.getHeight() > sourceTime)
				.filter(e -> e.hasSegments.cardinality() >= 1)				
				.collect(Collectors.toList());

		NetworkEdge iEdge = edges.get(Randomizer.nextInt(edges.size()));
		logHR -= Math.log(1.0/edges.size()); // Forward
		
		NetworkNode iParent = iEdge.parentNode;		
		NetworkNode iChild = iEdge.childNode;
		
        double delta = Math.abs(getDelta(iEdge.getLength()));        
        double logForwardDeltaProb = getLogDeltaProb(delta, iEdge.getLength());

        
        // get all potential reattachment Edges
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();

		// calculate the max height of the new child
		double minHeight = sourceTime;

		// check if delta is < iParent.getHeight()-sourceTime
		if (delta < (iParent.getHeight() - sourceTime)) {
			targetEdges.put(iEdge, sourceTime + delta);
		}else {
	    	targetEdges.putAll(getTargetEdgesUp(iParent.getParentEdges().get(0), delta, minHeight, iEdge));			
		}		
    	// get the other child edge that is not iEdge
    	logHR -= Math.log(1.0/targetEdges.size());
      				        	    	
        List<NetworkEdge> potentialNewTargets = targetEdges.keySet().stream().collect(Collectors.toList());        
		NetworkEdge jEdge = potentialNewTargets.get(Randomizer.nextInt(potentialNewTargets.size()));
				
		logHR += addReassortmentEdge(iEdge, sourceTime, jEdge, targetEdges.get(jEdge));
		
		if (logHR == Double.NEGATIVE_INFINITY) {
			return logHR; // no valid target edges
		}
	    	    
		int sizeAfter = (int) networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> e.hasSegments.cardinality() >= 1)				
				.count();

		logHR += Math.log(1.0/sizeAfter);  // Reverse
							    
		return logHR;
	}
	
	 
    private Map<NetworkEdge, Double> getTargetEdgesUp(NetworkEdge edge, double delta, double minHeight, NetworkEdge iEdge) {
    	
    	// initalize map
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();
        
        if (edge==iEdge)
        	return targetEdges; // do not move the edge itself
//        ii++;
  		if (edge.isRootEdge()) {
			double newHeight = edge.childNode.getHeight() + delta;
			targetEdges.put(edge, newHeight);			
			return targetEdges;
		}        
//		if (visitedNodes.contains(edge)) {
//			return targetEdges; // do not revisit nodes
//		}
//  		visitedNodes.add(edge);
  			
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
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(0), delta, minHeight, iEdge));
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(1), delta, minHeight, iEdge));				
			} else {
				// else it is a coalescent event
				targetEdges.putAll(getTargetEdgesUp(edge.parentNode.getParentEdges().get(0), delta, minHeight, iEdge));
				targetEdges.putAll(getTargetEdgesDown(getSisterEdge(edge), delta, minHeight, iEdge));
			}
		}
		return targetEdges;
	}
	
    private Map<NetworkEdge, Double> getTargetEdgesDown(NetworkEdge edge, double delta, double minHeight, NetworkEdge iEdge) {
        Map<NetworkEdge, Double> targetEdges = new HashMap<>();
		if (edge == iEdge)
			return targetEdges; // do not move the edge itself
        
        // get the difference between the current height and the new height
        delta -= edge.getLength();
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
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(0), delta, minHeight, iEdge));
			} else {
				// if coalescent event, follow both children down
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(0), delta, minHeight, iEdge));
				targetEdges.putAll(getTargetEdgesDown(edge.childNode.getChildEdges().get(1), delta, minHeight, iEdge));			
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
