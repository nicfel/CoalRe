package coalre.operators;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;


import beast.core.Input;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

public class NetworkExchange extends DivertSegmentOperator {
	final public Input<Boolean> isNarrowInput = new Input<>("isNarrow",
			"if true (default) a narrow exchange is performed, "
					+ "otherwise a wide exchange", true);


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
		return (n.getChildEdges().get(0).childNode.isCoalescence() ||
				n.getChildEdges().get(1).childNode.isCoalescence());
	}	


	/**
	 * Perform equivalent of narrow tree exchange on a network.
	 *
	 * @param	network
	 * @return	log of Hastings Ratio, or Double.NEGATIVE_INFINITY
	 * 			if proposal should not be accepted
	 */
	public double narrow(final Network network) {
		
		double logHR = 0.0;

		List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

		final List<NetworkEdge> possibleGrandParentEdges = networkEdges.stream()
				.filter(e -> e.childNode.isCoalescence())
				.filter(e -> hasCoalescenceKid(e.childNode))
				.collect(Collectors.toList());
		
		
		final int possibleGrandParents = possibleGrandParentEdges.size();
        if (possibleGrandParents < 1) {
            return Double.NEGATIVE_INFINITY;
        }
        logHR -= Math.log(1.0/possibleGrandParents);        

		final NetworkEdge grandParentEdge = possibleGrandParentEdges.
				get(Randomizer.nextInt(possibleGrandParents));
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

		if( !parent.isCoalescence()) {
			return Double.NEGATIVE_INFINITY;
		}

		
		final int childId = Randomizer.nextInt(parent.getChildEdges().size());
		final NetworkEdge childEdge = parent.getChildEdges().get(childId);

		
		logHR += exchangeEdges(childEdge, auntEdge, parent, grandParent);

		//TODO fix this workaround, too many states seem to be rejected  
		if (!allEdgesAncestral())
            return Double.NEGATIVE_INFINITY;
		
		networkEdges = new ArrayList<>(network.getEdges());
		
		final List<NetworkEdge> possibleGrandParentEdgesAfter = networkEdges.stream()
				.filter(e -> e.childNode.isCoalescence())
				.filter(e -> hasCoalescenceKid(e.childNode))
				.collect(Collectors.toList());
		
		final int possibleGrandParentsAfter = possibleGrandParentEdgesAfter.size();

		logHR += Math.log(1.0/possibleGrandParentsAfter);
		
		
		return logHR;
	}

	/**
	 * Perform equivalent of narrow tree exchange on a network.
	 *
	 * @param	network
	 * @return	log of Hastings Ratio, or Double.NEGATIVE_INFINITY
	 * 			if proposal should not be accepted
	 */
	public double wide(final Network network) {
		double logHR = 0.0;

		List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

		final List<NetworkEdge> possibleEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> !e.parentNode.isReassortment())
				.collect(Collectors.toList());
		
		final int nPossibleEdges = possibleEdges.size();
		logHR -= Math.log(1.0/(double)nPossibleEdges);

		final NetworkEdge iEdge = possibleEdges.
				get(Randomizer.nextInt(possibleEdges.size()));
		final NetworkNode i = iEdge.childNode;

		NetworkEdge jEdge = iEdge;

		while(jEdge == iEdge) {
			jEdge = possibleEdges.
					get(Randomizer.nextInt(possibleEdges.size()));
		}
		final NetworkNode j = jEdge.childNode;

		final NetworkNode p = iEdge.parentNode;
		final NetworkNode jP = jEdge.parentNode;


		if ((p != jP) && (i !=jP) && (j != p)
				&& (j.getHeight() < p.getHeight())
				&& (i.getHeight() < jP.getHeight())) {


			logHR += exchangeEdges(iEdge, jEdge, p, jP);
			
	        if (!allEdgesAncestral())
	            return Double.NEGATIVE_INFINITY;
			
			networkEdges = new ArrayList<>(network.getEdges());
			
			final List<NetworkEdge> possibleEdgesAfter = networkEdges.stream()
					.filter(e -> !e.isRootEdge())
					.filter(e -> !e.parentNode.isReassortment())
					.collect(Collectors.toList());
			
			final int nPossibleEdgesAfter = possibleEdgesAfter.size();
			logHR += Math.log(1.0/(double)nPossibleEdgesAfter);
	        
			return logHR;
		}
		else {
			return Double.NEGATIVE_INFINITY;
		}
	}


	/* exchange sub-nets whose root are i and j */
	protected double exchangeEdges(NetworkEdge iEdge, NetworkEdge jEdge,
			NetworkNode p, NetworkNode jP) {
		double logHR = 0.0;
		
		final NetworkEdge pEdge = p.getParentEdges().get(0);
		final NetworkEdge jPEdge = jP.getParentEdges().get(0);
		
		// After the exchange we want to add segments to the new ancestors
		// and remove from the old. Have to be careful not to remove segments
		// of siblings.

		final BitSet iSegs = iEdge.hasSegments;
		final BitSet jSegs = jEdge.hasSegments;

		final BitSet iSegsToRemove = (BitSet)iSegs.clone();
		iSegsToRemove.andNot(getSisterEdge(iEdge).hasSegments);
		iSegsToRemove.andNot(jSegs);

		final BitSet jSegsToRemove = (BitSet)jSegs.clone();
		jSegsToRemove.andNot(getSisterEdge(jEdge).hasSegments);
		jSegsToRemove.andNot(iSegs);
		
		final BitSet iSegsToAdd = (BitSet)iSegs.clone();
		iSegsToAdd.andNot(jSegs);
		
		final BitSet jSegsToAdd = (BitSet)jSegs.clone();
		jSegsToAdd.andNot(iSegs);
		
		p.removeChildEdge(iEdge);
		jP.removeChildEdge(jEdge);
		p.addChildEdge(jEdge);
		jP.addChildEdge(iEdge);
		
		logHR += removeSegmentsFromAncestors(jPEdge, jSegsToRemove);
		logHR += removeSegmentsFromAncestors(pEdge, iSegsToRemove);
		
		logHR -= addSegmentsToAncestors(jPEdge, iSegsToAdd);
		logHR -= addSegmentsToAncestors(pEdge, jSegsToAdd);
		
		return logHR;
	}


	/**
	 * Check that each edge is ancestral to at least one segment.
	 *
	 * @return true if all edges are ancestral.
	 */
	@Override
	public boolean allEdgesAncestral() {
		final Set<NetworkNode> nodeList = networkInput.get().getNodes();
		for (final NetworkNode node : nodeList) {
			for (final NetworkEdge parentEdge : node.getParentEdges()) {
				if (parentEdge.hasSegments.isEmpty())
					return false;
			}
		}
		return true;
	}

}
