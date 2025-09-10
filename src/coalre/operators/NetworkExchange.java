package coalre.operators;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import cern.colt.Arrays;
import coalre.distribution.NetworkEvent;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

public class NetworkExchange extends DivertSegmentOperator {
	final public Input<Boolean> isNarrowInput = new Input<>("isNarrow",
			"if true (default) a narrow exchange is performed, " + "otherwise a wide exchange", true);

	private boolean isNarrow;

	@Override
	public void initAndValidate() {
		super.initAndValidate();

		isNarrow = isNarrowInput.get();
	}

	@Override
	public double networkProposal() {

		double logHR;
		network.startEditing(this);

		if (isNarrow) {
			logHR = narrow(network);
		} else {
			logHR = wide(network);
		}
		return logHR;
	}

	private boolean hasCoalescenceKid(final NetworkNode n) {

		if (n.getChildEdges().get(0).childNode.isCoalescence() && n.getChildEdges().get(1).childNode.isCoalescence()) {
			return true;
		} else if (n.getChildEdges().get(0).childNode.isCoalescence()) {
			return n.getChildEdges().get(0).childNode.getHeight() > n.getChildEdges().get(1).childNode.getHeight();
		} else if (n.getChildEdges().get(1).childNode.isCoalescence()) {
			return n.getChildEdges().get(1).childNode.getHeight() > n.getChildEdges().get(0).childNode.getHeight();
		} else {
			return false;
		}
	}

	/**
	 * Perform equivalent of narrow tree exchange on a network.
	 *
	 * @param network
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should
	 *         not be accepted
	 */
	public double narrow(final Network network) {

		double logHR = 0.0;

//		List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

		final List<NetworkEdge> possibleGrandParentEdges = networkEdges.stream()
				.filter(e -> e.childNode.isCoalescence())
				.filter(e -> hasCoalescenceKid(e.childNode))
				.collect(Collectors.toList());

		final int possibleGrandParents = possibleGrandParentEdges.size();
		if (possibleGrandParents < 1) {
			return Double.NEGATIVE_INFINITY;
		}
		logHR -= Math.log(1.0 / possibleGrandParents);

		final NetworkEdge grandParentEdge = possibleGrandParentEdges.get(Randomizer.nextInt(possibleGrandParents));
		final NetworkNode grandParent = grandParentEdge.childNode;

		final List<NetworkEdge> possibleParentEdges = grandParent.getChildEdges();
		NetworkEdge parentEdge = possibleParentEdges.get(0);
		NetworkEdge auntEdge = possibleParentEdges.get(1);

		NetworkNode parent = parentEdge.childNode;
		NetworkNode aunt = auntEdge.childNode;

		if (parent.getHeight() < aunt.getHeight()) {
			auntEdge = possibleParentEdges.get(0);
			parentEdge = possibleParentEdges.get(1);

			parent = parentEdge.childNode;
			aunt = auntEdge.childNode;
		}

		if (!parent.isCoalescence()) {
			return Double.NEGATIVE_INFINITY;
		}

		final int childId = Randomizer.nextInt(parent.getChildEdges().size());
		final NetworkEdge childEdge = parent.getChildEdges().get(childId);

		logHR += exchangeEdges(childEdge, auntEdge, parent, grandParent);

		networkEdges = new ArrayList<>(network.getEdges());

		final List<NetworkEdge> possibleGrandParentEdgesAfter = networkEdges.stream()
				.filter(e -> e.childNode.isCoalescence()).filter(e -> hasCoalescenceKid(e.childNode))
				.collect(Collectors.toList());

		final int possibleGrandParentsAfter = possibleGrandParentEdgesAfter.size();

		logHR += Math.log(1.0 / possibleGrandParentsAfter);

		return logHR;
	}

	/**
	 * Perform equivalent of wide tree exchange on a network.
	 *
	 * @param network
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should
	 *         not be accepted
	 */
	public double wide(final Network network) {
		double logHR = 0.0;

		final List<NetworkEdge> possibleEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> e.hasSegments.cardinality() > 0)
				.filter(e -> isTrueCoalescence(e))
				.collect(Collectors.toList());

//		final int nPossibleEdges = possibleEdges.size();

		final NetworkEdge iEdge = possibleEdges.get(Randomizer.nextInt(possibleEdges.size()));

		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
		getTreeNodesDown(iEdge, (BitSet) iEdge.hasSegments.clone(), treeChildNodeList);

		// find all possible destination edges
		final List<NetworkEdge> possibleDestinationEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.getHeight()>iEdge.parentNode.getHeight())
				.filter(e -> e.childNode.getHeight()<iEdge.parentNode.getHeight())
				.filter(e -> segmentsOverlap(e.hasSegments, iEdge.hasSegments))
				.collect(Collectors.toList());
		
//		System.out.println(network);
//		System.out.println(iEdge.parentNode.getHeight());
		if (possibleDestinationEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;
		
		
		
		NetworkEdge	jEdge = possibleDestinationEdges.get(Randomizer.nextInt(possibleDestinationEdges.size()));
		
		logHR += addReassortmentEdge(iEdge, (iEdge.parentNode.getHeight()+iEdge.childNode.getHeight())/2, 
                jEdge, iEdge.parentNode.getHeight());
//		System.out.println(iEdge.parentNode.getParentEdges().get(0).hasSegments + " " + iEdge.parentNode.getParentEdges().get(1).hasSegments);
		if (iEdge.parentNode.getParentEdges().get(0).hasSegments.cardinality()==0)
			removeReassortmentEdge(iEdge.parentNode.getParentEdges().get(0));
		else
			removeReassortmentEdge(iEdge.parentNode.getParentEdges().get(1));
		
		int possibleEdgesRevers = (int) networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> e.hasSegments.cardinality() > 0)
				.filter(e -> isTrueCoalescence(e))
				.count();
		
		logHR += Math.log(possibleEdges.size()/(double)possibleEdgesRevers);

//		System.out.println(network);
//		System.exit(0);
		return logHR;
	}

	private boolean segmentsOverlap(BitSet hasSegments, BitSet hasSegments2) {
		for (int i = hasSegments.nextSetBit(0); i >= 0; i = hasSegments.nextSetBit(i + 1)) {
			if (hasSegments2.get(i))
				return true;
		}
		return false;
	}

	private boolean isTrueCoalescence(NetworkEdge e) {
		NetworkEdge sister = getSisterEdge(e);
		// loop over all segments in e, check if sister also has them, if so, return true
		for (int i = e.hasSegments.nextSetBit(0); i >= 0; i = e.hasSegments.nextSetBit(i + 1)) {
			if (sister.hasSegments.get(i))
				return true;
		}
		return false;
	}

	/* exchange sub-nets whose root are i and j */
	protected double exchangeEdges(NetworkEdge iEdge, NetworkEdge jEdge, NetworkNode p, NetworkNode jP) {
//		sanityCheck();

		double logHR = 0.0;

		final NetworkEdge pEdge = p.getParentEdges().get(0);
		final NetworkEdge jPEdge = jP.getParentEdges().get(0);

//		System.out.println(pEdge.getLength() + " " + jPEdge.getLength());
//		
//		System.out.println(network);

		// get all the segment trees involved
		int[] segmentIndices = iEdge.parentNode.segmentIndices;
		Integer[] treeNodesBelow = new Integer[network.getSegmentCount()];
		BitSet segsCoalBefore = (BitSet) getSisterEdge(iEdge).hasSegments.clone();
//		if (segmentIndices != null) {
		getTreeNodesDown(getSisterEdge(iEdge), segsCoalBefore, treeNodesBelow);
		segsCoalBefore.and(iEdge.hasSegments);
//		}

		// After the exchange we want to add segments to the new ancestors
		// and remove from the old. Have to be careful not to remove segments
		// of siblings.
		final BitSet iSegs = iEdge.hasSegments;
		final BitSet jSegs = jEdge.hasSegments;

		final BitSet iSegsToRemove = (BitSet) iSegs.clone();
		iSegsToRemove.andNot(getSisterEdge(iEdge).hasSegments);
		iSegsToRemove.andNot(jSegs);

		final BitSet jSegsToRemove = (BitSet) jSegs.clone();
		jSegsToRemove.andNot(getSisterEdge(jEdge).hasSegments);
		jSegsToRemove.andNot(iSegs);

		final BitSet iSegsToAdd = (BitSet) iSegs.clone();
		iSegsToAdd.andNot(jSegs);

		final BitSet jSegsToAdd = (BitSet) jSegs.clone();
		jSegsToAdd.andNot(iSegs);

		p.removeChildEdge(iEdge);
		jP.removeChildEdge(jEdge);
		p.addChildEdge(jEdge);
		jP.addChildEdge(iEdge);
//		System.out.println(network);
//		System.out.println(jEdge.parentNode.getParentEdges().get(0).getLength());
//		System.exit(0);

		NetworkEdge jEdgeParent = jEdge.parentNode.getParentEdges().get(0);

		jEdgeParent.hasSegments = (BitSet) jEdge.hasSegments.clone();
		jEdgeParent.hasSegments.or(getSisterEdge(jEdge).hasSegments);

//		logHR += removeSegmentsFromAncestors(jPEdge, jSegsToRemove);
//		logHR += removeSegmentsFromAncestors(pEdge, iSegsToRemove);
//		
//		logHR -= addSegmentsToAncestors(jPEdge, iSegsToAdd);
//		logHR -= addSegmentsToAncestors(pEdge, jSegsToAdd);

		// make a narrow move on all the segment indices
		BitSet segsCoalAfter = (BitSet) jEdge.hasSegments.clone();
		segsCoalAfter.and(getSisterEdge(jEdge).hasSegments);
//		System.out.println(network);

//		System.out.println(jEdge.parentNode.getHeight());

//		if (segmentIndices!=null) {

		BitSet segsToOperateOn = (BitSet) segsCoalBefore.clone();
		segsToOperateOn.or(segsCoalAfter);

//			System.out.println(segsCoalBefore + " " + segsCoalAfter + " " + segsToOperateOn);

		if (!segmentTrees.isEmpty()) {
			for (int i = 0; i < network.getSegmentCount(); i++) {
				if (!segsToOperateOn.get(i)) {
					continue; // segment not present in this tree
				}
				segmentTrees.get(i).startEditing(this);
	
				if (segsCoalBefore.get(i) && segsCoalAfter.get(i)) {
					Node parent = segmentTrees.get(i).getNode(segmentIndices[i]);
					Node grandParent = parent.getParent();
					Node uncle = getOtherChild(grandParent, parent);
					if (grandParent == null || parent == null) {
						throw new IllegalStateException("no root allowed: ");
					}
					// System.out.println(parent.isLeaf());
					// System.out.println(Arrays.toString(segmentIndices));
					Node child = getOtherChild(parent, segmentTrees.get(i).getNode(treeNodesBelow[i]));
	
					// System.out.println(segmentTrees.get(i) +";");
					exchangeNodes(child, uncle, parent, grandParent);
				} else if (segsCoalBefore.get(i)) {
					// resize the segment tree node
					int index = segmentIndices[i];
					Node nodeToReheight = segmentTrees.get(i).getNode(index);
					nodeToReheight.setHeight(jP.getHeight());
	
					if (jP.segmentIndices == null)
						jP.segmentIndices = new int[network.getSegmentCount()];
	
					jP.segmentIndices[i] = index;
					jEdge.parentNode.segmentIndices[i] = -1;
	
					// nodes coalesced before, but not after
	//					System.exit(0);
	//					return Double.NEGATIVE_INFINITY;
	
				} else {
					int index = jP.segmentIndices[i];
					Node nodeToReheight = segmentTrees.get(i).getNode(index);
					nodeToReheight.setHeight(jEdge.parentNode.getHeight());
	
					if (jEdge.parentNode.segmentIndices == null) {
						jEdge.parentNode.segmentIndices = new int[network.getSegmentCount()];
	//						System.out.println("init" + jEdge.parentNode.getHeight());
					}
	
					jEdge.parentNode.segmentIndices[i] = index;
					jP.segmentIndices[i] = -1;
	
	//					// nodes coalesced before, but not after
	//					System.out.println(network);
	//					System.out.println(i);
	//					System.out.println(jP.getHeight());
	//					System.out.println(jEdge.getLength());
	//					System.exit(0);
	
	//					return Double.NEGATIVE_INFINITY;
				}
	//		        System.out.println(segmentTrees.get(i).getLeafNodeCount());
	//		        System.out.println(segmentTrees.get(i) +";");
				// check if any edges have length <0
	//				for (Node n : segmentTrees.get(i).getExternalNodes()) {
	//					if (!n.isRoot()
	//							&& n.getLength() < 0) {
	//						System.err.println(i);
	//						throw new IllegalStateException("Negative height in segment tree: " + n.getHeight());
	//					}
	//				}
	//				for (Node n : segmentTrees.get(i).getInternalNodes()) {
	//					if (!n.isRoot()
	//							&& n.getLength() < 0) {
	//						System.err.println(i);
	//				        System.out.println(segmentTrees.get(i) +";");
	//						throw new IllegalStateException("Negative height in segment tree: " + n.getHeight());
	//					}
	//				}
	
	//				System.out.println(network);
	
	//		        System.exit(0);
	//		        
	//		        if (segmentTrees.get(i).getLeafNodeCount()!=350)
	//		        	System.exit(0);
			}
		}

//		}

//		System.out.println("exchanging edges: " + pEdge.childNode.getHeight() + " and " + jEdge.childNode.getHeight());
//		System.out.println(Arrays.toString(pEdge.parentNode.segmentIndices));
//		System.out.println(network);
//		System.exit(0);
//		sanityCheck();
		return logHR;
	}

	private void exchangeNodes(Node i, Node j, Node p, Node jP) {
		// precondition p -> i & jP -> j
		replace(p, i, j);
		replace(jP, j, i);
		// postcondition p -> j & p -> i
	}

	private void replace(final Node node, final Node child, final Node replacement) {
		node.removeChild(child);
		node.addChild(replacement);
		node.makeDirty(Tree.IS_FILTHY);
		replacement.makeDirty(Tree.IS_FILTHY);
	}

	private Node getOtherChild(final Node parent, final Node child) {
		if (parent.getLeft().getNr() == child.getNr()) {
			return parent.getRight();
		} else {
			return parent.getLeft();
		}
	}

	private boolean sanityCheck() {
		for (NetworkNode n : network.getNodes()) {
			if (n.isCoalescence()) {
				for (int i = 0; i < network.getSegmentCount(); i++) {
					if (n.getChildEdges().get(0).hasSegments.get(i) && n.getChildEdges().get(1).hasSegments.get(i)) {
						// check that the node height corresponds to the segment tree height
						double heightTree = 0;
						try {
							heightTree = segmentTrees.get(i).getNode(n.segmentIndices[i]).getHeight();
						} catch (Exception e) {
							System.out.println(n.getHeight());
							throw new IllegalStateException("Network node segment index not found in segment tree: "
									+ n.segmentIndices[i] + " for node " + n);
						}
						if (n.getHeight() != heightTree) {
							System.err.println("Network node height does not match segment tree height: "
									+ n.getHeight() + " != " + heightTree);
							System.err.println("Node: " + i);
							System.err.println("Segment tree: " + segmentTrees.get(i));
							throw new IllegalStateException("Network node height does not match segment tree height: "
									+ n.getHeight() + " != " + heightTree);
						}
					}
				}
			}
		}
		return true;
	}

	double addReassortmentEdge(NetworkEdge sourceEdge, double sourceTime, NetworkEdge destEdge, double destTime) {

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
		BitSet segsToDivert = (BitSet) sourceEdge.hasSegments.clone();
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


    


}
