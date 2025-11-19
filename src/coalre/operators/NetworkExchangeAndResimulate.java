package coalre.operators;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import coalre.distribution.NetworkEvent;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

public class NetworkExchangeAndResimulate extends SubNetworkLeapAndResimulate {
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

		if (isNarrow) {
			logHR = narrow(network);
		} else {
			logHR = wide_old(network);
		}
		
		return logHR;
	}

	private int isg(final Node n) {
		if (n.getLeft().isLeaf() && n.getRight().isLeaf()) {
			return 0;
		}
		
		if (n.getLeft().getHeight() < n.getRight().getHeight()) {
			if (n.getRight().isLeaf()) {
				return 0;
			}
		}else {
			if (n.getLeft().isLeaf()) {
				return 0;
			}
		}
		
		
		return 1;
	}

	private int sisg(final Node n) {
		return n.isLeaf() ? 0 : isg(n);
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
		// pick a random segment
		int segment = Randomizer.nextInt(network.getSegmentCount());
		// get a node on the tree that is eligble for a narrow exchange
		final int internalNodes = segmentTrees.get(segment).getInternalNodeCount();
		if (internalNodes <= 1) {
			return Double.NEGATIVE_INFINITY;
		}
//		System.out.println(segmentTrees.get(segment) +";");
		
		int validGP = 0;
		{
			for (int i = internalNodes + 1; i < 1 + 2 * internalNodes; ++i) {
				validGP += isg(segmentTrees.get(segment).getNode(i));
			}
		}
		
		if (validGP == 0) {
			return Double.NEGATIVE_INFINITY;
		}

		Node grandParent = segmentTrees.get(segment).getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
		while (isg(grandParent)==0) {
			grandParent = segmentTrees.get(segment).getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
		}

		Node parentIndex = grandParent.getLeft();
		Node uncle = grandParent.getRight();
		if (parentIndex.getHeight() < uncle.getHeight()) {
			parentIndex = grandParent.getRight();
			uncle = grandParent.getLeft();
		}

		if (parentIndex.isLeaf()) {
			throw new IllegalStateException("Unexpected leaf node selected in");
		}


//		final int c2 = sisg(parentIndex) + sisg(uncle);

		int parentInd = parentIndex.getNr();
		int uncleInd = uncle.getNr();
		// find the network edge that has parentIndex coalescing and uncle coalescing
		List<NetworkEdge> parentEdge = networkEdges.stream().filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence()).filter(e -> e.parentNode.segmentIndices != null)
				.filter(e -> e.parentNode.segmentIndices[segment] == parentInd)
				.filter(e -> e.parentNode.getChildEdges().get(0).hasSegments.get(segment))
				.filter(e -> e.parentNode.getChildEdges().get(1).hasSegments.get(segment))
				.collect(Collectors.toList());

		List<NetworkEdge> targetEdge = networkEdges.stream().filter(e -> e.childNode.segmentIndices != null)
				.filter(e -> e.childNode.segmentIndices[segment] == uncleInd)
				.filter(e -> e.childNode.isLeaf()
						|| (e.childNode.getChildEdges().get(0).hasSegments.get(segment)
								&& e.childNode.getChildEdges().get(1).hasSegments.get(segment)))
				.collect(Collectors.toList());

		if (parentEdge.size() != 2) {
			return Double.NEGATIVE_INFINITY;
		}
		NetworkEdge iEdge = parentEdge.get(Randomizer.nextInt(2));


		// find the edge above targetEdge, following segment, that has parent>iEdge.parent
		// and child<iEdge.parent
		NetworkEdge jEdge = getEdgeWithHeight(targetEdge.get(0), iEdge.parentNode.getHeight(), segment);
		
				
		BitSet segsToDivert = (BitSet) iEdge.hasSegments.clone();
		
		logHR += divertSegmentsWithResimulation(iEdge, jEdge, segsToDivert, iEdge.parentNode.getHeight());

		
		int validGPafter = 0;
		{
			for (int i = internalNodes + 1; i < 1 + 2 * internalNodes; ++i) {
				validGPafter += isg(segmentTrees.get(segment).getNode(i));
			}
		}

//		final int validGPafter = validGP - c2 + sisg(parentIndex) + sisg(uncle);

		logHR += Math.log((float) validGP / validGPafter);
//		System.out.println(segmentTrees.get(segment));
//		System.out.println(logHR);
//		System.out.println("validGP: "+validGP);
//		System.out.println("validGPafter: "+validGPafter);
//		System.out.println(Math.log((float) validGP / validGPafter));
//
//		if (logHR>5) {
//			System.exit(0);
//		}
		return logHR;
	}

	private NetworkEdge getEdgeWithHeight(NetworkEdge edge, double height, int segment) {
		if (edge.parentNode.getHeight() > height && edge.childNode.getHeight() < height) {
			return edge;
		} else {
			for (NetworkEdge parentEdge : edge.parentNode.getParentEdges()) {
				if (parentEdge.hasSegments.get(segment)) {
					return getEdgeWithHeight(parentEdge, height, segment);
				}
			}
		}
		return null;
	}

	
	int reaDiff = 0;
	int count =0;
	double logHRTotal =0.0;
	
	/**
	 * Perform equivalent of wide tree exchange on a network.
	 *
	 * @param network
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should
	 *         not be accepted
	 */
	public double wide_old(final Network network) {
		double logHR = 0.0;

		final List<NetworkEdge> possibleEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
//				.filter(e -> isCoalNode(e))
				.collect(Collectors.toList());

		int edgeCount = possibleEdges.size();

//		reaDiff = (int) networkEdges.stream()
//				.filter(e -> !e.isRootEdge())
//				.filter(e -> e.parentNode.isReassortment())
//				.count();

		final NetworkEdge iEdge = possibleEdges.get(Randomizer.nextInt(possibleEdges.size()));

		// find all possible destination edges
		final List<NetworkEdge> possibleDestinationEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.getHeight() > iEdge.parentNode.getHeight())
				.filter(e -> e.childNode.getHeight() < iEdge.parentNode.getHeight())
//				.filter(e -> segmentOverlap(e, iEdge.hasSegments))
				.collect(Collectors.toList());

		if (possibleDestinationEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;

//		int forwardDestCount = possibleDestinationEdges.size();

//		System.out.println(network.getExtendedNewick(0));


		NetworkEdge jEdge = possibleDestinationEdges.get(Randomizer.nextInt(possibleDestinationEdges.size()));
		BitSet segsToDivert = (BitSet) iEdge.hasSegments.clone();
		logHR += divertSegmentsWithResimulation(iEdge, jEdge, segsToDivert, iEdge.parentNode.getHeight());
//		System.out.println(logHR);


		final List<NetworkEdge> possibleEdgesReverse = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
//				.filter(e -> isCoalNode(e))
				.collect(Collectors.toList());

//		reaDiff -= (int) networkEdges.stream()
//				.filter(e -> !e.isRootEdge())
//				.filter(e -> e.parentNode.isReassortment())
//				.count();


		int reverseEdgeCount = possibleEdgesReverse.size();
		logHR += Math.log((double) edgeCount / reverseEdgeCount);

		// Account for destination edge selection
		// In reverse, we would select a destination edge at the height of jEdge.parentNode
		// Note: After the move, jEdge now has the segments that were on iEdge
//		final List<NetworkEdge> possibleReverseDestinationEdges = networkEdges.stream()
//				.filter(e -> !e.isRootEdge())
//				.filter(e -> e.parentNode.getHeight() > jEdge.parentNode.getHeight())
//				.filter(e -> e.childNode.getHeight() < jEdge.parentNode.getHeight())
//				.filter(e -> segmentOverlap(e, segsToDivert))
//				.collect(Collectors.toList());

//		if (possibleReverseDestinationEdges.isEmpty())
//			return Double.NEGATIVE_INFINITY;
//
//		int reverseDestCount = possibleReverseDestinationEdges.size();
//		logHR += Math.log((double) forwardDestCount / reverseDestCount);

//		if (logHR == Double.POSITIVE_INFINITY) {
//			System.out.println(network.getExtendedNewick(0));
//			System.exit(0);
//		}

//		logHRTotal += logHR;
//		count++;
//		System.out.println("Reassortment difference: " + ((double) reaDiff/count) + " logHR: " + (logHRTotal/count));
		return logHR;
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
	
	private boolean segmentOverlap(NetworkEdge e, BitSet segs) {
		for (int i = segs.nextSetBit(0); i >= 0; i = segs.nextSetBit(i + 1)) {
			if (e.hasSegments.get(i)) {
				return true;
			}
		}
		return false;
	}



//	public double wide(final Network network) {
//		double logHR = 0.0;
//		// pick a random segment
//		int segment = Randomizer.nextInt(network.getSegmentCount());
//		// get a node on the tree that is eligble for a narrow exchange
//		final int internalNodes = segmentTrees.get(segment).getInternalNodeCount();
//		if (internalNodes <= 1) {
//			return Double.NEGATIVE_INFINITY;
//		}
//		
//		Node grandParent = segmentTrees.get(segment).getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
//		while (grandParent.isRoot()) {
//			grandParent = segmentTrees.get(segment).getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
//		}
//
//		Node parentIndex = grandParent.getLeft();
//		if (Randomizer.nextBoolean()) {
//			parentIndex = grandParent.getRight();
//		}
//
////		final int c2 = sisg(parentIndex) + sisg(uncle);
//
//		int parentInd = parentIndex.getNr();
//		// find the network edge that has parentIndex coalescing and uncle coalescing
//		List<NetworkEdge> parentEdge = networkEdges.stream().filter(e -> !e.isRootEdge())
//				.filter(e -> e.parentNode.isCoalescence()).filter(e -> e.parentNode.segmentIndices != null)
//				.filter(e -> e.parentNode.segmentIndices[segment] == parentInd)
//				.filter(e -> e.parentNode.getChildEdges().get(0).hasSegments.get(segment))
//				.filter(e -> e.parentNode.getChildEdges().get(1).hasSegments.get(segment))
//				.collect(Collectors.toList());
//
//		// find all co-existing edges that have the segment
//
//		if (parentEdge.size() != 2) {
//			return Double.NEGATIVE_INFINITY;
//		}
//		NetworkEdge iEdge = parentEdge.get(Randomizer.nextInt(2));
//
//
//		// find the edge above targetEdge, following segment, that has parent>iEdge.parent
//		// and child<iEdge.parent
//		NetworkEdge jEdge = getEdgeWithHeight(targetEdge.get(0), iEdge.parentNode.getHeight(), segment);
//		
//				
//		BitSet segsToDivert = (BitSet) iEdge.hasSegments.clone();
//		
//		logHR += divertSegmentsWithResimulation(iEdge, jEdge, segsToDivert, iEdge.parentNode.getHeight());
//
//		
//		int validGPafter = 0;
//		{
//			for (int i = internalNodes + 1; i < 1 + 2 * internalNodes; ++i) {
//				validGPafter += isg(segmentTrees.get(segment).getNode(i));
//			}
//		}
//
////		final int validGPafter = validGP - c2 + sisg(parentIndex) + sisg(uncle);
//
//		logHR += Math.log((float) validGP / validGPafter);
////		System.out.println(segmentTrees.get(segment));
////		System.out.println(logHR);
////		System.out.println("validGP: "+validGP);
////		System.out.println("validGPafter: "+validGPafter);
////		System.out.println(Math.log((float) validGP / validGPafter));
////
////		if (logHR>5) {
////			System.exit(0);
////		}
//		return logHR;
//	}



}
