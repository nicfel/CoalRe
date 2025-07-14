package coalre.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import cern.colt.Arrays;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class DivertSegmentOperator extends EmptyEdgesNetworkOperator {

	public Input<Boolean> divertOneSegmentInput = new Input<>("divertOneSegment",
			"If true, only one segment is diverted", false);

	@Override
	public double networkProposal() {
		double logHR = 0.0;

		List<NetworkEdge> sourceEdges = network.getEdges().stream().filter(e -> e.childNode.isReassortment())
				.filter(e -> e.hasSegments.cardinality() > 0).collect(Collectors.toList());

		if (sourceEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;

		logHR -= Math.log(1.0 / sourceEdges.size());

		NetworkEdge sourceEdge = sourceEdges.get(Randomizer.nextInt(sourceEdges.size()));
		NetworkEdge destEdge = getSpouseEdge(sourceEdge);

		BitSet segsToDivert;

		if (divertOneSegmentInput.get()) {
			segsToDivert = getRandomSegment(sourceEdge.hasSegments);
			logHR -= Math.log(1. / sourceEdge.hasSegments.cardinality());
		} else {
			segsToDivert = getRandomUnconditionedSubset(sourceEdge.hasSegments);
			logHR -= getLogUnconditionedSubsetProb(sourceEdge.hasSegments);
		}

		// set all
		for (int i = 0; i < segmentTrees.size(); i++) {
//			if (!segsToDivert.get(i)) {
				segmentsChanged.set(i, false);
//			} else {
//				segmentsChanged.set(i, true);
//			}
		}

		if (segsToDivert.cardinality() == 0)
			return Double.NEGATIVE_INFINITY;

		network.startEditing(this);

//        Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
//        System.out.println(segsToDivert);
//        getTreeNodesDown(sourceEdge, segsToDivert, treeChildNodeList);
//        System.out.println(segsToDivert);
		// the start nodes tell us the original parent, the other child, find the new
		// target node and height
		// update the segments
		logHR += divertSegments(destEdge, sourceEdge, segsToDivert);

//        List<NetworkNode> targetNodeIndices = new ArrayList<>();
//        for (int i =0; i < network.getSegmentCount();i++)
//        	targetNodeIndices.add(null);
//        Integer[] newTargetNode = new Integer[network.getSegmentCount()];
//        Double[] newNodeHeights = new Double[network.getSegmentCount()];
//
//        getTreeNodes(destEdge, segsToDivert, targetNodeIndices, newTargetNode, newNodeHeights);

//        System.out.println(Arrays.toString(treeChildNodeList));
//        System.out.println(Arrays.toString(newTargetNode));
//        System.out.println(Arrays.toString(newNodeHeights));
//        System.out.println(targetNodeIndices);
//        
//        // update the corresponding tree
//        for (int i =0; i < network.getSegmentCount();i++) {
//        	if (treeChildNodeList[i]!=null) {
//        		System.out.println(targetNodeIndices.get(i).getHeight());
//        		
//        	}
//    	}        
//        targetNodeIndices.clear();
//        System.exit(0);

		if (divertOneSegmentInput.get()) {
			logHR += Math.log(1.0 / destEdge.hasSegments.cardinality());
		} else {
			logHR += getLogUnconditionedSubsetProb(destEdge.hasSegments);
		}

		int reverseSourceEdgeCount = (int) (network.getEdges().stream().filter(e -> e.childNode.isReassortment())
				.filter(e -> e.hasSegments.cardinality() > 0).count());

		logHR += Math.log(1.0 / reverseSourceEdgeCount);

		return logHR;
	}

	protected double divertSegments(NetworkEdge destEdge, NetworkEdge sourceEdge, BitSet segsToDivert) {
		double logHR = 0.0;
		System.out.println(segsToDivert + " " + sourceEdge.childNode.getHeight());
		System.out.println(network);
		
		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
		
		
		getTreeNodesDown(sourceEdge, segsToDivert, treeChildNodeList);
//		System.out.println(segsToDivert);

		logHR -= addSegmentsToAncestors(destEdge, segsToDivert);
		logHR += removeSegmentsFromAncestors(sourceEdge, segsToDivert);

		List<NetworkNode> targetNodeIndices = new ArrayList<>();
		for (int i = 0; i < network.getSegmentCount(); i++)
			targetNodeIndices.add(null);
		Integer[] newTargetNode = new Integer[network.getSegmentCount()];
		Double[] newNodeHeights = new Double[network.getSegmentCount()];

		getTreeNodes(destEdge, segsToDivert, targetNodeIndices, newTargetNode, newNodeHeights);

//      System.out.println(Arrays.toString(treeChildNodeList));
//      System.out.println(Arrays.toString(newTargetNode));
//      System.out.println(Arrays.toString(newNodeHeights));
//      System.out.println(targetNodeIndices);
//      
		System.out.println(network);

		// update the corresponding tree
		for (int i = 0; i < network.getSegmentCount(); i++) {
			if (treeChildNodeList[i] != null) {
//				System.out.println(i+ " " + targetNodeIndices.get(i).getHeight());
//				System.out.println(segmentTrees.get(i) +";");
				Node startChild = segmentTrees.get(i).getNode(treeChildNodeList[i]);
				Node oldParent = startChild.getParent();
				Node newSibling = segmentTrees.get(i).getNode(newTargetNode[i]);

				if (oldParent.isRoot()) {

					if (oldParent==newSibling) {
						System.out.println("lsssala");

						oldParent.setHeight(newNodeHeights[i]);
					}else {
						// new root
						throw new IllegalArgumentException("lala");
					}
					
				}else {
					Node oldGrandParent = oldParent.getParent();
					Node oldSibling = oldParent.getChild(0).getNr() == treeChildNodeList[i] ? oldParent.getChild(1) : oldParent.getChild(0);
					
//					System.out.println(oldGrandParent.getHeight());
					oldGrandParent.removeChild(oldParent);
					oldGrandParent.addChild(oldSibling);
					Node newGrandParent = newSibling.getParent();
					if (newGrandParent==null) {
						oldParent.setParent(null);
//						System.out.println(newSibling.isRoot());
						oldParent.removeChild(oldSibling);
						oldParent.addChild(newSibling);
						oldParent.setHeight(newNodeHeights[i]);
//						System.out.println(oldParent.getParent().getHeight());
						segmentTrees.get(i).setRoot(oldParent);
						System.out.println("lala");

					}else {
						newGrandParent.removeChild(newSibling);
						oldParent.removeChild(oldSibling);
						newGrandParent.addChild(oldParent);
						oldParent.addChild(newSibling);
						oldParent.setHeight(newNodeHeights[i]);
					}
				}

				
				
				
				System.out.println(segmentTrees.get(i) +";");
//				System.exit(0);
			}
		}

		return logHR;
	}

	private BitSet getRandomSegment(BitSet hasSegments) {
		// pick a random segment to divert.
		BitSet segsToDivert = new BitSet();
		segsToDivert.clear();

		int segToDivert = Randomizer.nextInt(hasSegments.cardinality());
		int segIdx = hasSegments.nextSetBit(0);
		for (int i = 0; i < segToDivert; i++) {
			segIdx = hasSegments.nextSetBit(segIdx + 1);
		}
		segsToDivert.set(segIdx);
		return segsToDivert;
	}

//    double removeSegmentsFromAncestors(NetworkEdge edge, BitSet segsToRemove) {
//    	int[] treeNodesIn = new int[network.getSegmentCount()];
//    	Node[] treeNodes = new Node[network.getSegmentCount()];
//    	
//    	return removeSegmentsFromAncestors(edge, segsToRemove, treeNodes, treeNodesIn);
//    }

	void getTreeNodesDown(NetworkEdge edge, BitSet segsToRemove, Integer[] treeNodeList) {
		segsToRemove = (BitSet) segsToRemove.clone();
		segsToRemove.and(edge.hasSegments);

		if (segsToRemove.isEmpty())
			return;
		
		if (edge.childNode.isReassortment()) {
			getTreeNodesDown(edge.childNode.getChildEdges().get(0), segsToRemove, treeNodeList);
		} else if (edge.childNode.isCoalescence()) { // check if this node corresponds to the tree node
			if (segmentTrees.size() > 0 && edge.childNode.segmentIndices != null) {
				for (int segIdx = 0; segIdx < segsToRemove.length(); segIdx++) {
					// get the corresponding tree node
					if (segsToRemove.get(segIdx) && edge.childNode.segmentIndices[segIdx] != -1) {
						treeNodeList[segIdx] = edge.childNode.segmentIndices[segIdx];
						// set bitset to false
						segsToRemove.set(segIdx, false);
					}
				}
			}

			getTreeNodesDown(edge.childNode.getChildEdges().get(0), segsToRemove, treeNodeList);
			getTreeNodesDown(edge.childNode.getChildEdges().get(1), segsToRemove, treeNodeList);

		} else {
			// loop over the segsToRemove Left and update the Tree Nodes
			for (int segIdx = 0; segIdx < segsToRemove.length(); segIdx++) {
				if (segsToRemove.get(segIdx)) {
					treeNodeList[segIdx] = edge.childNode.segmentIndices[segIdx];
				}
			}
		}

		return;
	}

	void getTreeNodes(NetworkEdge edge, BitSet segsToRemove, List<NetworkNode> targetNodeIndices,
			Integer[] treeNodeList, Double[] nodeHeightList) {
		segsToRemove = (BitSet) segsToRemove.clone();
		segsToRemove.and(edge.hasSegments);

		if (segsToRemove.isEmpty())
			return;

		if (edge.isRootEdge())
			return;

		if (edge.parentNode.isReassortment()) {

			getTreeNodes(edge.parentNode.getParentEdges().get(0), segsToRemove, targetNodeIndices, treeNodeList,
					nodeHeightList);
			getTreeNodes(edge.parentNode.getParentEdges().get(1), segsToRemove, targetNodeIndices, treeNodeList,
					nodeHeightList);

		} else {
			if (segmentTrees.size() > 0) {
				BitSet sibSegs = getSisterEdge(edge).hasSegments;

				if (segmentTrees.size() > 0) {
					for (int segIdx = 0; segIdx < segsToRemove.length(); segIdx++) {
						if (segsToRemove.get(segIdx) && sibSegs.get(segIdx)) {
							nodeHeightList[segIdx] = edge.parentNode.getHeight();
							treeNodeList[segIdx] = getTreeNodeIndex(getSisterEdge(edge), segIdx);
							targetNodeIndices.set(segIdx, edge.parentNode);
							segsToRemove.set(segIdx, false);
						}
					}
				}

			}

			getTreeNodes(edge.parentNode.getParentEdges().get(0), segsToRemove, targetNodeIndices, treeNodeList,
					nodeHeightList);

		}

		return;
	}

	private Integer getTreeNodeIndex(NetworkEdge edge, int segIdx) {
		if (edge.childNode.isLeaf()) {
			return edge.childNode.segmentIndices[segIdx];
		} else if (edge.childNode.isReassortment()) {
			return getTreeNodeIndex(edge.childNode.getChildEdges().get(0), segIdx);
		} else {
			if (edge.childNode.segmentIndices != null && edge.childNode.segmentIndices[segIdx] != -1) {
				return edge.childNode.segmentIndices[segIdx];
			}
			return edge.childNode.getChildEdges().get(0).hasSegments.get(segIdx)
					? getTreeNodeIndex(edge.childNode.getChildEdges().get(0), segIdx)
					: getTreeNodeIndex(edge.childNode.getChildEdges().get(1), segIdx);
		}
	}

	/**
	 * Remove segments from this edge and ancestors.
	 *
	 * @param edge         edge at which to start removal
	 * @param segsToRemove segments to remove from edge and ancestors
	 * @return log probability of reverse operation
	 */
	private double removeSegmentsFromAncestors(NetworkEdge edge, BitSet segsToRemove) {
		double logP = 0.0;

		segsToRemove = (BitSet) segsToRemove.clone();
		segsToRemove.and(edge.hasSegments);

		if (segsToRemove.isEmpty())
			return logP;

		edge.hasSegments.andNot(segsToRemove);

		if (edge.isRootEdge())
			return logP;

		if (edge.parentNode.isReassortment()) {

			logP += Math.log(0.5) * segsToRemove.cardinality();

			logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove);
			logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(1), segsToRemove);

		} else {
//        	if (segmentTrees.size() > 0) {
//	        	BitSet segsToRemoveLeft = (BitSet) segsToRemove.clone();
//	        	segsToRemoveLeft.and(getSisterEdge(edge).hasSegments);
//	        	// loop over the segsToRemove Left and update the Tree Nodes
//	        	while (segsToRemoveLeft.nextSetBit(0) != -1) {
//					int segIdx = segsToRemoveLeft.nextSetBit(0);
//					// get the corresponding tree node
//					Node treeNode = segmentTrees.get(segIdx).getNode(edge.parentNode.segmentIndices[segIdx]);
//					treeNodes[segIdx] = treeNode;
//					// get the child node index that is does not have number treeNodesIn[segIdx]
//					Node otherChild = treeNode.getChild(0).getNr;
//							
//				}	
//        	}
			segsToRemove.andNot(getSisterEdge(edge).hasSegments);

			logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove);

		}

		return logP;
	}

	/**
	 * Add segments to this edge and ancestors.
	 *
	 * @param edge      edge at which to start addition
	 * @param segsToAdd segments to add to the edge and ancestors
	 * @return log probability of operation
	 */
	private double addSegmentsToAncestors(NetworkEdge edge, BitSet segsToAdd) {
		double logP = 0.0;

		segsToAdd = (BitSet) segsToAdd.clone();
		segsToAdd.andNot(edge.hasSegments);

		if (segsToAdd.isEmpty())
			return logP;

		edge.hasSegments.or(segsToAdd);

		if (edge.isRootEdge())
			return logP;

		if (edge.parentNode.isReassortment()) {

			BitSet segsToAddLeft = new BitSet();
			BitSet segsToAddRight = new BitSet();

			for (int segIdx = segsToAdd.nextSetBit(0); segIdx != -1; segIdx = segsToAdd.nextSetBit(segIdx + 1)) {
				if (Randomizer.nextBoolean())
					segsToAddLeft.set(segIdx);
				else
					segsToAddRight.set(segIdx);

				logP += Math.log(0.5);
			}

			logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAddLeft);
			logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(1), segsToAddRight);

		} else {

			logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAdd);
		}

		return logP;
	}

	protected BitSet getRandomUnconditionedSubset(BitSet sourceSegments) {
		BitSet destSegments = new BitSet();

		destSegments.clear();

		for (int segIdx = sourceSegments.nextSetBit(0); segIdx != -1; segIdx = sourceSegments.nextSetBit(segIdx + 1)) {

			if (Randomizer.nextBoolean())
				destSegments.set(segIdx);
		}

		return destSegments;
	}

	protected double getLogUnconditionedSubsetProb(BitSet sourceSegments) {
		return sourceSegments.cardinality() * Math.log(0.5);
	}

}
