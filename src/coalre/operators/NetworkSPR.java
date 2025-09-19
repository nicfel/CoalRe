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

public class NetworkSPR extends DivertSegmentOperator {
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
		if (divertOneSegmentInput.get())
			return oneseg(network);
		else
			return general(network);
	}

	/**
	 * Perform equivalent of wide tree exchange on a network.
	 *
	 * @param network
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should
	 *         not be accepted
	 */
	public double general(final Network network) {
		double logHR = 0.0;

		final List<NetworkEdge> possibleEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> e.childNode.isReassortment())
				.filter(e -> !e.parentNode.getParentEdges().get(0).isRootEdge())
				.collect(Collectors.toList());

		if (possibleEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;

		final NetworkEdge iEdge = possibleEdges.get(Randomizer.nextInt(possibleEdges.size()));

		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
		getTreeNodesDown(iEdge, (BitSet) iEdge.hasSegments.clone(), treeChildNodeList);

		// find all possible destination edges
		final List<NetworkEdge> possibleDestinationEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.getHeight()>iEdge.childNode.getHeight())
				.collect(Collectors.toList());
		
		if (possibleDestinationEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;
		
		NetworkEdge	jEdge = possibleDestinationEdges.get(Randomizer.nextInt(possibleDestinationEdges.size()));
		
		double oldLower = Math.max(iEdge.childNode.getHeight(), getSisterEdge(iEdge).childNode.getHeight());
		double oldDiff = iEdge.parentNode.getParentEdges().get(0).parentNode.getHeight() - oldLower;
		
		// randomly sample a time
		double lower = Math.max(jEdge.childNode.getHeight(), iEdge.childNode.getHeight());
		double diff = jEdge.parentNode.getHeight() - lower;
		double attatchmentTime = Randomizer.nextDouble() * diff + lower;
		logHR += Math.log(diff); 
		logHR -= Math.log(oldDiff);
				
		logHR += addReassortmentEdge(iEdge, (iEdge.parentNode.getHeight()+iEdge.childNode.getHeight())/2, 
                jEdge, attatchmentTime);
		
		if (iEdge.parentNode.getParentEdges().get(0).hasSegments.cardinality()==0)
			removeReassortmentEdge(iEdge.parentNode.getParentEdges().get(0));
		else
			removeReassortmentEdge(iEdge.parentNode.getParentEdges().get(1));
		
		
		return logHR;
	}

	
	public double oneseg(final Network network) {
		double logHR = 0.0;
		
		List<Integer> possibleEdges = new ArrayList<>();
		for (int i = 0; i < networkEdges.size(); i++) {
			NetworkEdge e = networkEdges.get(i);
			if (!e.isRootEdge() 
					&& e.parentNode.isCoalescence() 
					&& e.hasSegments.cardinality() == 1
					&& !e.parentNode.getParentEdges().get(0).isRootEdge()) {
				int segIdx = e.hasSegments.nextSetBit(0);
				if (getSisterEdge(e).hasSegments.get(segIdx))
						possibleEdges.add(i);
			}
		}
		if (possibleEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;
			

		final NetworkEdge iEdge = networkEdges.get(possibleEdges.get(Randomizer.nextInt(possibleEdges.size())));
		int segment = iEdge.hasSegments.nextSetBit(0);
		

		Integer[] treeChildNodeList = new Integer[network.getSegmentCount()];
		getTreeNodesDown(iEdge, (BitSet) iEdge.hasSegments.clone(), treeChildNodeList);

		// find all possible destination edges
		final List<NetworkEdge> possibleDestinationEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.get(segment))
				.filter(e -> e.parentNode.getHeight()>iEdge.childNode.getHeight())
				.collect(Collectors.toList());
		
		if (possibleDestinationEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;
		
		NetworkEdge	jEdge = possibleDestinationEdges.get(Randomizer.nextInt(possibleDestinationEdges.size()));
		
		double oldLower = Math.max(iEdge.childNode.getHeight(), getSisterEdge(iEdge).childNode.getHeight());
		double oldDiff = iEdge.parentNode.getParentEdges().get(0).parentNode.getHeight() - oldLower;
		
		// randomly sample a time
		double lower = Math.max(jEdge.childNode.getHeight(), iEdge.childNode.getHeight());
		double diff = jEdge.parentNode.getHeight() - lower;
		double attatchmentTime = Randomizer.nextDouble() * diff + lower;
		logHR += Math.log(diff); 
		logHR -= Math.log(oldDiff);
				
		logHR += addReassortmentEdge(iEdge, (iEdge.parentNode.getHeight()+iEdge.childNode.getHeight())/2, 
                jEdge, attatchmentTime);
		
		if (iEdge.parentNode.getParentEdges().get(0).hasSegments.cardinality()==0)
			removeReassortmentEdge(iEdge.parentNode.getParentEdges().get(0));
		else
			removeReassortmentEdge(iEdge.parentNode.getParentEdges().get(1));
		
		
		return logHR;
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
