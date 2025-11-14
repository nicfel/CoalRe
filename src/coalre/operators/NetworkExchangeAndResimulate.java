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
		network.startEditing(this);

		for (int i = 0; i < segmentTrees.size(); i++) {
			segmentsChanged.set(i, false);
		}

		if (isNarrow) {
			logHR = narrow(network);
		} else {
			logHR = wide(network);
		}
		
		return logHR;
	}

	private int isg(final Node n) {
		return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
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

		Node grandParent = segmentTrees.get(segment).getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
		while (grandParent.getLeft().isLeaf() && grandParent.getRight().isLeaf()) {
			grandParent = segmentTrees.get(segment).getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
		}

		Node parentIndex = grandParent.getLeft();
		Node uncle = grandParent.getRight();
		if (parentIndex.getHeight() < uncle.getHeight()) {
			parentIndex = grandParent.getRight();
			uncle = grandParent.getLeft();
		}

		if (parentIndex.isLeaf()) {
			// tree with dated tips
			return Double.NEGATIVE_INFINITY;
		}

		int validGP = 0;
		{
			for (int i = internalNodes + 1; i < 1 + 2 * internalNodes; ++i) {
				validGP += isg(segmentTrees.get(segment).getNode(i));
			}
		}

		final int c2 = sisg(parentIndex) + sisg(uncle);

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

		
		
		
		// Get network events to determine lineages at different times
		List<NetworkEvent> networkEventList = coalescentDistr.intervals.getNetworkEventList();
				
		BitSet segsToDivert = (BitSet) iEdge.hasSegments.clone();
		
		logHR += divertSegmentsWithResimulation(iEdge, jEdge, segsToDivert, iEdge.parentNode.getHeight());

		
		
		final int validGPafter = validGP - c2 + sisg(parentIndex) + sisg(uncle);

		logHR += Math.log((float) validGP / validGPafter);
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

	/**
	 * Perform equivalent of wide tree exchange on a network.
	 *
	 * @param network
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should
	 *         not be accepted
	 */
	public double wide(final Network network) {
		double logHR = 0.0;

		final List<NetworkEdge> possibleEdges = networkEdges.stream().filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence()).collect(Collectors.toList());
		
		int edgeCount = possibleEdges.size();

		final NetworkEdge iEdge = possibleEdges.get(Randomizer.nextInt(possibleEdges.size()));

		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
		getTreeNodesDown(iEdge, (BitSet) iEdge.hasSegments.clone(), treeChildNodeList);

		// find all possible destination edges
		final List<NetworkEdge> possibleDestinationEdges = networkEdges.stream().filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.getHeight() > iEdge.parentNode.getHeight())
				.filter(e -> e.childNode.getHeight() < iEdge.parentNode.getHeight()).collect(Collectors.toList());

		if (possibleDestinationEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;

		NetworkEdge jEdge = possibleDestinationEdges.get(Randomizer.nextInt(possibleDestinationEdges.size()));
		BitSet segsToDivert = (BitSet) iEdge.hasSegments.clone();
		logHR += divertSegmentsWithResimulation(iEdge, jEdge, segsToDivert, iEdge.parentNode.getHeight());
		
		

		final List<NetworkEdge> possibleEdgesReverse = networkEdges.stream().filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence()).collect(Collectors.toList());
		

		
		int reverseEdgeCount = possibleEdgesReverse.size();
		logHR += Math.log((double) edgeCount / reverseEdgeCount);
		return logHR;
	}


}
