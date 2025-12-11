package coalre.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import coalre.distribution.CoalescentWithReassortment;
import coalre.distribution.NetworkEvent;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.statistics.NetworkStatsLogger;

import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

public class BranchUniform extends DivertSegmentOperator {
	
    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<Boolean> gaussianInput = new Input<>("gaussian", "Gaussian (=true=default) or uniform delta", true);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, " +
            "i.e. -1. (when positive, gets multiplied by network-height/log2(n-taxa).", -1.0);


    CoalescentWithReassortment coalescentDistr;

    boolean useMaxHeight = false;
    
    double addProb = 0.5;
    
    double size;
    double limit;

    
    @Override
    public void initAndValidate() {
        super.initAndValidate();

		size = sizeInput.get();
		limit = limitInput.get();
    }

    @Override
    public double networkProposal() {
        double logHR = 0.0;
        
        List<Integer> candidateEdges = new ArrayList<>();
		for (int i = 0; i < networkEdges.size(); i++) {
			NetworkEdge edge = networkEdges.get(i);
			if (!edge.isRootEdge() && edge.childNode.isReassortment()
					&& edge.parentNode.isCoalescence() 
					&& !edge.parentNode.getParentEdges().get(0).isRootEdge()
					&& edge.hasSegments.cardinality() == 1
					&& validEdge(edge)) {
				candidateEdges.add(i);
			}
		}
//		System.out.println("candidate edges: " + candidateEdges.size());
		if (candidateEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;
		
		
//        System.out.println(network);
		

		
        NetworkEdge edgeToMove = networkEdges.get(candidateEdges.get(Randomizer.nextInt(candidateEdges.size())));
        int segment = edgeToMove.hasSegments.nextSetBit(0);
        logHR -= Math.log(1.0/((double) candidateEdges.size()));
        
        // for this edge, get the corresponding tree node
        Node n = segmentTrees.get(segment).getNode(edgeToMove.parentNode.segmentIndices[segment]);
        // get the minimum and maximum times we can move this edge to
        double minTime = Math.max(n.getLeft().getHeight(),n.getRight().getHeight());
        double maxTime = n.isRoot() ? network.getRootEdge().childNode.getHeight() : n.getParent().getHeight();
        
        // randomly sample the new desttim
        double newDestTime = Randomizer.nextDouble() * (maxTime - minTime) + minTime;

        // find out if the path of the reassortment edge corresponds to left or right, by making a new method that tracks the 
        // edge to move path down until segment index matches the edge to move
        int index = getChildNodeIndex(edgeToMove, segment);
        Node correspondingChild = n.getLeft().getNr() == index ? n.getLeft() : n.getRight();


        NetworkEdge oldDestEdge = getSisterEdge(edgeToMove);
        NetworkEdge oldSourceEdge = edgeToMove.childNode.getChildEdges().get(0);
        
        double oldSourceTime = edgeToMove.childNode.getHeight();
        double oldDestTime = edgeToMove.parentNode.getHeight();
//        System.out.println(segmentTrees.get(segment) +";");        

		removeReassortmentEdge(edgeToMove);
        // check in again on the corresponding child node and sample the new time
		double sourceTimeMax = Math.min(correspondingChild.getParent().getHeight(), newDestTime);
		
        double newSourceTime = Randomizer.nextDouble() * (sourceTimeMax - correspondingChild.getHeight()) + correspondingChild.getHeight();

        logHR -= Math.log(1.0/(sourceTimeMax - correspondingChild.getHeight()));
        
		double reversSourceTimeMax = Math.min(correspondingChild.getParent().getHeight(), oldDestTime);
        logHR += Math.log(1.0/(reversSourceTimeMax - correspondingChild.getHeight()));
//        System.out.println(newSourceTime + " " + oldSourceTime);
//        System.out.println(newDestTime + " " + oldDestTime);

//		System.out.println(oldDestEdge.isRootEdge());
		// System.out.println(network);
		// System.exit(0);
		
		
//		System.out.println(network);
        
    	List<NetworkEdge> targetEdgesDest = new ArrayList<>();
    	List<NetworkEdge> targetEdgesSource = new ArrayList<>();

        if (newDestTime>oldDestTime)
       	    logHR += getTargetEdgesUp(oldDestEdge, newDestTime, targetEdgesDest, segment);
        else
            logHR += getTargetEdgesDown(oldDestEdge, newDestTime, targetEdgesDest, segment);

        if (newSourceTime>oldSourceTime) 
            logHR += getTargetEdgesUp(oldSourceEdge, newSourceTime, targetEdgesSource, segment);
        else
            logHR += getTargetEdgesDown(oldSourceEdge, newSourceTime, targetEdgesSource, segment);

//    	if (delta>0) {
//        	// follow the segments up
//        	logHR += getTargetEdgesUp(oldDestEdge, newDestTime, targetEdgesDest, segment);
//        	logHR += getTargetEdgesUp(oldSourceEdge, newSourceTime, targetEdgesSource, segment);
//        }else {
//        	// follow the segments down
//        	logHR += getTargetEdgesDown(oldDestEdge, newDestTime, targetEdgesDest, segment);
//        	logHR += getTargetEdgesDown(oldSourceEdge, newSourceTime, targetEdgesSource, segment);
//        }
//        
    	if (targetEdgesDest.isEmpty() || targetEdgesSource.isEmpty())
    		return Double.NEGATIVE_INFINITY;
    	
        // pick one of the target edges at random
        NetworkEdge newDestEdge = targetEdgesDest.get(Randomizer.nextInt(targetEdgesDest.size()));
        NetworkEdge newSourceEdge = targetEdgesSource.get(Randomizer.nextInt(targetEdgesSource.size()));
        
        logHR -= Math.log(1.0/((double) targetEdgesDest.size()));
        logHR -= Math.log(1.0/((double) targetEdgesSource.size()));
        
        // make the calculation for the reverse HR
        List<NetworkEdge> reverseTargetEdgesDest = new ArrayList<>();
        List<NetworkEdge> reverseTargetEdgesSource = new ArrayList<>();

        if (newDestTime>oldDestTime)
            logHR += getTargetEdgesDown(newDestEdge, oldDestTime, reverseTargetEdgesDest, segment);
        else
            logHR += getTargetEdgesUp(newDestEdge, oldDestTime, reverseTargetEdgesDest, segment);

        if (newSourceTime>oldSourceTime)
            logHR += getTargetEdgesUp(newSourceEdge, oldSourceTime, reverseTargetEdgesSource, segment);
        else
            logHR += getTargetEdgesDown(newSourceEdge, oldSourceTime, reverseTargetEdgesSource, segment);



//        if (delta>0) {
//        	logHR += getTargetEdgesDown(newDestEdge, oldDestTime, reverseTargetEdgesDest, segment);
//        	logHR += getTargetEdgesDown(newSourceEdge, oldSourceTime, reverseTargetEdgesSource, segment);
//        }else {
//			logHR += getTargetEdgesUp(newDestEdge, oldDestTime, reverseTargetEdgesDest, segment);
//			logHR += getTargetEdgesUp(newSourceEdge, oldSourceTime, reverseTargetEdgesSource, segment);
//        }
        
        logHR += Math.log(1.0/((double) reverseTargetEdgesDest.size()));
        logHR += Math.log(1.0/((double) reverseTargetEdgesSource.size()));
//        System.out.println(network);
        
        addReassortmentEdge(newSourceEdge, newSourceTime, newDestEdge, newDestTime, segment);
      
        List<Integer> reverseCandidateEdges = new ArrayList<>();
		for (int i = 0; i < networkEdges.size(); i++) {
			NetworkEdge edge = networkEdges.get(i);
			if (!edge.isRootEdge()&& edge.childNode.isReassortment()
					&& edge.parentNode.isCoalescence() 
					&& !edge.parentNode.getParentEdges().get(0).isRootEdge()
					&& edge.hasSegments.cardinality() == 1
					&& validEdge(edge)) {
				reverseCandidateEdges.add(i);
			}
		}
		
        logHR += Math.log(1.0/((double) reverseCandidateEdges.size()));
        
        // check if any edges have length<=0 if so print         	System.out.println(network);

        for (NetworkEdge e: networkEdges) {
        	if (!e.isRootEdge() &&e.getLength()<=0) {
				System.err.println("edge length <=0");
//				System.err.println(saveOldNetwork);
				System.err.println(network);
				System.exit(0);
			}
    	}
        
        // do reverse delta calculation
//        double probabReverse = getLogDeltaProb(-delta, currentDistance + delta);
//        logHR += probabReverse - probabForward;
        
        
//        System.out.println(segmentTrees.get(segment) +";\n");        

//        System.out.println(network+"\n");
        
        return logHR;
    }
    boolean validEdge(NetworkEdge e) {
    	return getSisterEdge(e).hasSegments.get(e.hasSegments.nextSetBit(0));
    }

    int getChildNodeIndex(NetworkEdge edge, int segment) {
        if (edge.childNode.isLeaf()) {
            return edge.childNode.segmentIndices[segment];
        }else if (edge.childNode.isCoalescence()) {
            if (edge.childNode.getChildEdges().get(0).hasSegments.get(segment)
                && edge.childNode.getChildEdges().get(1).hasSegments.get(segment)) {
                    return edge.childNode.segmentIndices[segment];
            }else{
                if (edge.childNode.getChildEdges().get(0).hasSegments.get(segment))
                    return getChildNodeIndex(edge.childNode.getChildEdges().get(0), segment);
                else
                    return getChildNodeIndex(edge.childNode.getChildEdges().get(1), segment);
            }
        } else {
            return getChildNodeIndex(edge.childNode.getChildEdges().get(0), segment);
        }
    }

    double addReassortmentEdge(NetworkEdge sourceEdge, double sourceTime,
                               NetworkEdge destEdge, double destTime, int segment) {

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
    	BitSet segsToDivert = new BitSet();
    	segsToDivert.set(segment);
            
        logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
        networkEdges.add(reassortmentEdge);
        networkEdges.add(newEdge1);
        networkEdges.add(newEdge2);


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
    
    
    // Traverse up the network (forward in time) - only direct ancestors
    private double getTargetEdgesUp(NetworkEdge currentEdge, double destTime, 
                                List<NetworkEdge> targetEdges, int segment) {
        
    	double logHR = 0.0;
    	
        if (currentEdge.isRootEdge()) {
            // Reached root
            targetEdges.add(currentEdge);
            return logHR;
        }
        
		if (currentEdge.childNode.getHeight() < destTime && currentEdge.parentNode.getHeight() >= destTime) {
			// Current edge spans destTime - add to target edges
			targetEdges.add(currentEdge);
			return logHR;
		}else {
            for (NetworkEdge e : currentEdge.parentNode.getParentEdges()) {
            	if (e.hasSegments.get(segment))
            		logHR += getTargetEdgesUp(e, destTime, targetEdges, segment);
            }

		}
        return logHR;
    }

	private double getTargetEdgesDown(NetworkEdge currentEdge, double destTime, 
                                List<NetworkEdge> targetEdges, int segment) {
        
    	double logHR = 0.0;
    	
		if (currentEdge.isRootEdge() && currentEdge.childNode.getHeight() <= destTime) {
			// Reached root
			targetEdges.add(currentEdge);
			return logHR;
		}else if (currentEdge.childNode.getHeight() <= destTime && currentEdge.parentNode.getHeight() > destTime) {
    		// Current edge spans destTime - add to target edges
            targetEdges.add(currentEdge);
            return logHR;
    	}
    	
        if (currentEdge.isLeafEdge()) {
            return logHR;
        }
        
        for (NetworkEdge e : currentEdge.childNode.getChildEdges()) {
        	if (e.hasSegments.get(segment))
        		logHR += getTargetEdgesDown(e, destTime, targetEdges, segment);
        }
        return logHR;
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