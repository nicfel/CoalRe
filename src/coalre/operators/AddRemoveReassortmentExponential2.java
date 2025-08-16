package coalre.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import coalre.distribution.CoalescentWithReassortment;
import coalre.distribution.NetworkEvent;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.statistics.NetworkStatsLogger;
import test.beastfx.app.inputeditor.BaseTest.e;

import java.util.*;
import java.util.stream.Collectors;

public class AddRemoveReassortmentExponential2 extends DivertSegmentOperator {

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

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
        size = sizeInput.get();
        
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
			try {
				logHR = removeRecombination();
			} catch (Exception e) {
				return Double.NEGATIVE_INFINITY;
			}
		}		
		return logHR;
	}

	double addRecombination() {
			
		double logHR = 0.0;
		int segment = Randomizer.nextInt(network.getSegmentCount());
//		if (divertOneSegmentInput.get()) {
	    logHR -= Math.log(1.0 / network.getSegmentCount()); // pick random segment
//		}


		Tree segmentTree = segmentTrees.get(segment);

		// Get the segment tree for the chosen segment
		int sourceNodeIndex = Randomizer.nextInt(segmentTree.getNodeCount());
		while (segmentTree.getNode(sourceNodeIndex).isRoot()) 
			sourceNodeIndex = Randomizer.nextInt(segmentTree.getNodeCount());
				
		logHR -= Math.log(1.0 / (segmentTree.getNodeCount()-1)); // pick random source node
		
		Node iNode = segmentTree.getNode(sourceNodeIndex);
		Node otherChild = iNode.getParent().getLeft()==iNode ? iNode.getParent().getRight() : iNode.getParent().getLeft();
		
		
        double delta = Math.abs(getDelta(iNode.getLength()));        
        logHR -= getLogDeltaProb(delta, iNode.getLength()); //pick delta

        // get all potential reattachment Edges
        Map<Node, Double> targetEdges = new HashMap<>();
		// calculate the max height of the new child
		double minHeight = iNode.getHeight();

    	targetEdges.putAll(getTargetEdgesUp(iNode.getParent(), delta, minHeight));
    	// get the other child edge that is not iEdge
    	targetEdges.putAll(getTargetEdgesDown(otherChild, delta, minHeight));
    	
    	logHR -= Math.log(1.0/targetEdges.size());

        List<Node> potentialNewTargets = targetEdges.keySet().stream().collect(Collectors.toList());        
        Node destNode = potentialNewTargets.get(Randomizer.nextInt(potentialNewTargets.size()));
        double destTime = targetEdges.get(destNode);
        double timeDiff = Math.min(iNode.getParent().getHeight(), destTime) - iNode.getHeight();
        
        double sourceTime = Randomizer.nextDouble() * timeDiff + iNode.getHeight();
        logHR -= Math.log(1.0/timeDiff);
        
        
        // find the source edge and the dest edge   
        Optional<NetworkEdge> sourceEdgeStart = findNetworkEdge(sourceNodeIndex, segment, iNode.getHeight());
        Optional<NetworkEdge> destEdgeStart = findNetworkEdge(destNode.getNr(), segment, destNode.getHeight());
        
        // check that both edges have been found, otherwise, throw an exception
		if (!sourceEdgeStart.isPresent() || !destEdgeStart.isPresent()) {
			throw new IllegalArgumentException("Source or destination edge not found for recombination.");
		}
		
		// find the correct source and destEdge as the the edges themeselves or their parent for which 
		// child node < time and parent node > time
		NetworkEdge sourceEdge = getParentEdge(sourceEdgeStart.get(), sourceTime, segment);
		NetworkEdge destEdge = getParentEdge(destEdgeStart.get(), destTime, segment);
		
		// Create new reassortment edge
		logHR += addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime, segment);
		// check if one of the parents of source edge hs cardinality 0
//		if (sourceEdge.parentNode.getParentEdges().get(0).hasSegments.cardinality() == 0) {
//			System.out.println(network);
//			removeReassortmentEdge(sourceEdge.parentNode.getParentEdges().get(0), segment);
//		} else if(sourceEdge.parentNode.getParentEdges().get(1).hasSegments.cardinality() == 0) {
//			removeReassortmentEdge(sourceEdge.parentNode.getParentEdges().get(1), segment);
//		}
		

		if (logHR == Double.NEGATIVE_INFINITY)
			return Double.NEGATIVE_INFINITY;

		if (divertOneSegmentInput.get()) {
			// HR contribution for reverse move
			int nRemovableEdges = (int) networkEdges.stream()
					.filter(e -> !e.isRootEdge())
					.filter(e -> e.childNode.isReassortment())
					.filter(e -> e.parentNode.isCoalescence())
					.filter(e -> e.hasSegments.cardinality() == 1)
					.filter(e -> actualCoalescent(e))
					.count();
	
			logHR += Math.log(1.0 / nRemovableEdges);
		}else {
			
			int nRemovableEdges = (int) networkEdges.stream()
					.filter(e -> !e.isRootEdge())
					.filter(e -> e.childNode.isReassortment())
					.filter(e -> e.parentNode.isCoalescence())
					.filter(e -> isCoalNode(e))
					.count();
			logHR += Math.log(1.0/nRemovableEdges);
			
			// calculate the coalSegment
			BitSet coalSegments = (BitSet) destEdge.hasSegments.clone();
			coalSegments.and(getSisterEdge(destEdge).hasSegments);
			
			
			logHR += Math.log(1.0/coalSegments.cardinality());
			

		}
		return logHR;
	}



	double removeRecombination() {
	    double logHR = 0.0;

	    // Find removable edges (same as before)
	    List<Integer> removableEdges = new ArrayList<>();
	    int segment =-1;
	    NetworkEdge edgeToRemove;
	    if (divertOneSegmentInput.get()) {
        	    
		    for (int i = 0; i < networkEdges.size(); i++) {
		        NetworkEdge edge = networkEdges.get(i);
		        if (!edge.isRootEdge() 
		        		&& edge.childNode.isReassortment() 
		        		&& edge.parentNode.isCoalescence()
		                && edge.hasSegments.cardinality() == 1
		                && actualCoalescent(edge)) {
		            removableEdges.add(i);
		        }
		    }
		    
		    // Sample edge to remove
		    edgeToRemove = networkEdges.get(removableEdges.get(Randomizer.nextInt(removableEdges.size())));
		    segment = edgeToRemove.hasSegments.nextSetBit(0);

	    }else {

		    for (int i = 0; i < networkEdges.size(); i++) {
		        NetworkEdge edge = networkEdges.get(i);
		        if (!edge.isRootEdge() 
		        		&& edge.childNode.isReassortment() 
		        		&& edge.parentNode.isCoalescence()
		                && isCoalNode(edge)) {
		            removableEdges.add(i);
		        }
		    }
		    
		    edgeToRemove = networkEdges.get(removableEdges.get(Randomizer.nextInt(removableEdges.size())));
		    // get all edges that coalesce at this node
		    BitSet coalSegments = (BitSet) edgeToRemove.hasSegments.clone();
		    coalSegments.and(getSisterEdge(edgeToRemove).hasSegments);
		    		    
		    segment = coalSegments.nextSetBit(Randomizer.nextInt(coalSegments.cardinality()));
		    
	    	// The HR contribution should account for choosing this segment from ALL segments in forward move
	    	logHR -= Math.log(1.0 / coalSegments.cardinality()); // pick random segment
	    }


	    if (removableEdges.isEmpty()) {
	        return Double.NEGATIVE_INFINITY;
	    }
	    
	    logHR -= Math.log(1.0 / removableEdges.size());


	    // Extract parameters from the current configuration
	    double oldDestTime = edgeToRemove.parentNode.getHeight();
	    
	    NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
	    NetworkEdge destEdge = getSisterEdge(edgeToRemove);

	    if (destEdge.childNode == edgeToRemove.childNode)
	        destEdge = sourceEdge;
	    
	    
	    // Find the source and target nodes in the segment tree
	    NetworkEdge sourceStartEdge = getSegmentNodeIndex(segment, sourceEdge);
	    int sourceNodeIndex = sourceStartEdge.childNode.segmentIndices[segment];
	    
	    // Find the target node by following the dest edge
	    NetworkEdge destStartEdge = getSegmentNodeIndex(segment, destEdge);
	    int destNodeIndex = destStartEdge.childNode.segmentIndices[segment];
	    
	    // get the common ancestor time of the two tree nodes
	    Tree segmentTree = segmentTrees.get(segment);
	    Node startNode = segmentTree.getNode(sourceNodeIndex);
	    Node destNodeBefore = segmentTree.getNode(destNodeIndex);
	    
	    double prevRootHeight = segmentTree.getRoot().getHeight();
	    boolean wasRoot = startNode.getParent().isRoot();
	    
//    	System.out.println(segmentTree +";");
	    logHR += removeReassortmentEdge(edgeToRemove, segment);
//	    System.out.println(network);
//	    System.exit(0);

	    if (logHR == Double.NEGATIVE_INFINITY) {
	        return Double.NEGATIVE_INFINITY;
	    }

	    Node destNodeAfter = startNode.getParent().getLeft()==startNode ? startNode.getParent().getRight() : startNode.getParent().getLeft();
	    double delta = 0.0;
	    // if start node or destNodebefore are root, do something else
	    if (wasRoot) {
	    	delta = startNode.getParent().getHeight() - prevRootHeight;
	    } else if (startNode.getParent().isRoot()) {
			// is newly created root
	    	delta = oldDestTime - startNode.getParent().getHeight();
		}else {
		    double mrcaTime = getMRCAtime(startNode, destNodeBefore);
		    double pathLength = 2 * mrcaTime - startNode.getHeight() - oldDestTime;
		    delta = pathLength - startNode.getLength();
		}
	    
	    // compute the implied delta
	    delta = Math.abs(delta);
	    
	    logHR += Math.log(1.0 / (segmentTree.getNodeCount() - 1));
		
		logHR += getLogDeltaProb(delta, startNode.getLength());
		
		Node otherChild = startNode.getParent().getLeft() == startNode ? startNode.getParent().getRight() : startNode.getParent().getLeft();
		double minHeight = startNode.getHeight();

		Map<Node, Double> targetEdges = new HashMap<>();
    	targetEdges.putAll(getTargetEdgesUp(startNode.getParent(), delta, minHeight));
    	// get the other child edge that is not iEdge
    	targetEdges.putAll(getTargetEdgesDown(otherChild, delta, minHeight));
    	
    	logHR += Math.log(1.0 / targetEdges.size());
    	
    	// contribution for the source time
        double timeDiff = Math.min(startNode.getParent().getHeight(), oldDestTime) - startNode.getHeight();
        logHR += Math.log(1.0 / timeDiff);
    	
//        if (divertOneSegmentInput.get()) {
            logHR += Math.log(1.0 / network.getSegmentCount());
//        }


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

		if (divertOneSegmentInput.get()) {
			BitSet segsToDivert = new BitSet();
			segsToDivert.clear();
			segsToDivert.set(segment);
			logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
		}else {
			BitSet segsToDivert = new BitSet();
			segsToDivert.clear();
			segsToDivert.set(segment);	
			
			// for each other segment, add it with 0.5 prob
			for (int i = 0; i < network.getSegmentCount();i++) {
				if (i != segment && newEdge1.hasSegments.get(i)) {
					if (Randomizer.nextBoolean()) {
						segsToDivert.set(i);
					}
					logHR -= Math.log(0.5);

				}
			}
			logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
		}

		networkEdges.add(reassortmentEdge);
		networkEdges.add(newEdge1);
		networkEdges.add(newEdge2);
		return logHR;
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
		
        if (!divertOneSegmentInput.get()) {
        	// In removeReassortmentEdge(), replace lines 439-448 with:
        	NetworkEdge sourceEdgeInReverse = edgeToRemove.childNode.getChildEdges().get(0);
        	int availableForRandomSelection = 0;

        	for (int i = 0; i < network.getSegmentCount(); i++) {
        	    if (i != segment && sourceEdgeInReverse.hasSegments.get(i)) {
        	        availableForRandomSelection++;
                	logHR += Math.log(0.5);
        	    }
        	}
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

    
    private Map<Node, Double> getTargetEdgesUp(Node node, double delta, double minHeight) {
    	
    	// initalize map
        Map<Node, Double> targetEdges = new HashMap<>();
  		if (node.isRoot()) {
			double newHeight = node.getHeight() + delta;
			targetEdges.put(node, newHeight);			
			return targetEdges;
		}        
  			
        // get the difference between the current height and the new height
        delta -= node.getLength();

		if (delta <= 0.0) {
			double newHeight = node.getHeight() + delta + node.getLength();

			targetEdges.put(node, newHeight);			
			return targetEdges; 
		} else {
			// else it is a coalescent event
			targetEdges.putAll(getTargetEdgesUp(node.getParent(), delta, minHeight));
			Node otherChild = node.getParent().getLeft() == node ? node.getParent().getRight() : node.getParent().getLeft();
			targetEdges.putAll(getTargetEdgesDown(otherChild, delta, minHeight));
		}
		return targetEdges;
	}

	
    private Map<Node, Double> getTargetEdgesDown(Node node, double delta, double minHeight) {
        Map<Node, Double> targetEdges = new HashMap<>();

        // get the difference between the current height and the new height
        delta -= node.getLength();

        if (delta <= 0.0) {
        	double newHeight = node.getHeight() - delta - node.getLength();
            if (newHeight < minHeight)
            	return targetEdges; // no valid targets
            
			targetEdges.put(node, node.getParent().getHeight() - (delta + node.getLength()));
			return targetEdges; 
		} else {
			if (node.getHeight() < minHeight) {
				return targetEdges; // no valid targets
			}

			// proceed to the next event
			if (node.isLeaf()) 
				return targetEdges; // no valid targets
			
			// if coalescent event, follow both children down
			targetEdges.putAll(getTargetEdgesDown(node.getLeft(), delta, minHeight));
			targetEdges.putAll(getTargetEdgesDown(node.getRight(), delta, minHeight));			
			
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

}
