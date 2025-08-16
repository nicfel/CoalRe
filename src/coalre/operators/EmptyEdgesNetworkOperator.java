package coalre.operators;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

public abstract class EmptyEdgesNetworkOperator extends NetworkOperator {

	public Input<Double> emptyAlphaInput = new Input<>("emptyAlpha",
			"Mean of exponential used for choosing root attachment times.", 0.1);

	public Input<Double> lambdaInput = new Input<>("lambda",
			"lambda of the poisson distribution for how many empty edges to add.", 0.1);

	public Input<Boolean> addRemoveEmptyEdgesInput = new Input<>("addRemoveEmptyEdges",
			"adds empty edges before calling the networkproposal and then removes all empty edges at the end again",
			true);
	
	public Input<Integer> removedEdgesOffset = new Input<>("removedEdgesOffset",
			"Offset for the removed edges. This is used to avoid the operator to remove edges that were added in the same proposal.",
			0);
		
	private double emptyAlpha;
	private double lambda;
	protected List<NetworkEdge> networkEdges;
	int netChange;
	int ii=0;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		emptyAlpha = emptyAlphaInput.get();
		lambda = lambdaInput.get();
		segmentsChanged = new BitSet(segmentTrees.size());
	}

	@Override
	public double proposal() {

		double logHR = 0.0;

		network.startEditing(this);

		networkEdges = new ArrayList<>(network.getEdges());
//		
//		
//		
//		netChange = 0;
//		double logHRproposal=0.0;
//		// Adds empty network edges
//		if (addRemoveEmptyEdgesInput.get()) {
//			logHRproposal = addEmptyNetworkSegments();
//			logHR += logHRproposal;
//
//			if (logHR == Double.NEGATIVE_INFINITY)
//				return Double.NEGATIVE_INFINITY;
//		}

       	logHR += networkProposal();
//
//
//		if (addRemoveEmptyEdgesInput.get()) {
//
//			if (logHR == Double.NEGATIVE_INFINITY)
//				return Double.NEGATIVE_INFINITY;
//			try {
//				
//				// remove empty reassortment edges
//				logHR += RemoveAllEmptyNetworkSegments();
//			} catch (Exception e) {
//				// if there are no empty reassortment edges, this can happen when removing empty
//				// edges. This was previously taken care of
//				// by the all edges ancestral check, but that took too much time
//				return Double.NEGATIVE_INFINITY;
//			}
//		}
//
//		if (logHR == Double.POSITIVE_INFINITY)
//			return Double.NEGATIVE_INFINITY;
//
//		// check that the networkEdges size is the same as newly calculating them
////		if (networkEdges.size()!=network.getEdges().size()) {
////			System.err.println(networkEdges.size() + " != " + network.getEdges().size() + " after proposal" +  this.getClass());
////			throw new IllegalArgumentException("Network edges size changed during proposal, this should not happen! "
////                    + "Please report this as a bug to the developers of coalre.");
////		}
//		
////		if (!networkTerminatesAtMRCA())
////			return Double.NEGATIVE_INFINITY;
//		
//		ii++;
//		if (netChange <-5) {
//			System.out.println(this.getID() + " " + netChange + " " + logHR + " " + logHRproposal);
//		}

		return logHR;
	}

	private double addEmptyNetworkSegments() {
		double logHR = 0.0;

		// randomly sample the number of edges to add
		int nrEmptyEdges = (int) Randomizer.nextPoisson(lambda);
		netChange+= nrEmptyEdges;

		for (int i = 0; i < nrEmptyEdges; i++) {
			logHR += addEmptyReassortment();
		}
//    	if (this.getID().equals("Co")) {
//    		System.out.println("added " + nrEmptyEdges + " empty edges");
//    	}

		logHR -= Math.log(Math.pow(lambda, nrEmptyEdges)) - lambda - logFactorial(nrEmptyEdges);

		return logHR;
	}

	double addEmptyReassortment() {
		double logHR = 0.0;

		// add empty reassortment edges to non empty edges
//		
//		List<NetworkEdge> possibleSourceEdges = networkEdges.stream().filter(e -> !e.isRootEdge())
//				.collect(Collectors.toList());
				

		NetworkEdge sourceEdge = networkEdges.get(Randomizer.nextInt(networkEdges.size()));
		while (sourceEdge.isRootEdge())
			sourceEdge = networkEdges.get(Randomizer.nextInt(networkEdges.size()));
		
		double sourceTime = Randomizer.nextDouble() * sourceEdge.getLength() + sourceEdge.childNode.getHeight();

		logHR -= Math.log(1.0 / (double) (networkEdges.size()-1)) + Math.log(1.0 / sourceEdge.getLength());

		List<Integer> possibleDestEdges = new ArrayList<>();
		for (int i = 0; i < networkEdges.size(); i++) {
			if (!networkEdges.get(i).isRootEdge() && networkEdges.get(i).parentNode.getHeight() >= sourceTime) {
				possibleDestEdges.add(i);
			}
			if (networkEdges.get(i).isRootEdge()) {
				possibleDestEdges.add(i);
			}
		}

		NetworkEdge destEdge = networkEdges.get(possibleDestEdges.get(Randomizer.nextInt(possibleDestEdges.size())));

		// works
		logHR -= Math.log(1.0 / possibleDestEdges.size());

		double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

		double destTime;
		if (destEdge.isRootEdge()) {
			// works
			destTime = minDestTime + Randomizer.nextExponential(1.0 / emptyAlpha);
			logHR -= -(1.0 / emptyAlpha) * (destTime - minDestTime) + Math.log(1.0 / emptyAlpha);
		} else {
			destTime = Randomizer.nextDouble() * (destEdge.parentNode.getHeight() - minDestTime) + minDestTime;
			logHR -= Math.log(1.0 / (destEdge.parentNode.getHeight() - minDestTime));
		}

		// Create new reassortment edge
		logHR += addEmptyReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);

		if (logHR == Double.NEGATIVE_INFINITY)
			return Double.NEGATIVE_INFINITY;

		// HR contribution for reverse move
		int nRemovableEdges = 0;
		for (NetworkEdge e : networkEdges) {
			if (!e.isRootEdge() && e.childNode.isReassortment()
					&& e.parentNode.isCoalescence()  && e.hasSegments.cardinality() == 0) {
				nRemovableEdges++;
			}
		}
		
		// works
		logHR += Math.log(1.0 / nRemovableEdges);

		return logHR;
	}

	// Only adds the reassortment edge, but does not diverge segments
	double addEmptyReassortmentEdge(NetworkEdge sourceEdge, double sourceTime, NetworkEdge destEdge, double destTime) {

		double logHR = 0.0;

		NetworkNode sourceNode = new NetworkNode();
		sourceNode.setHeight(sourceTime);

		NetworkNode oldSourceEdgeParent = sourceEdge.parentNode;
		oldSourceEdgeParent.removeChildEdge(sourceEdge);
		sourceNode.addChildEdge(sourceEdge);

		NetworkEdge newEdge1 = new NetworkEdge();
		sourceNode.addParentEdge(newEdge1);
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

		// update edges
		networkEdges.add(newEdge1);
		networkEdges.add(newEdge2);
		networkEdges.add(reassortmentEdge);

		return logHR;
	}

	public double RemoveAllEmptyNetworkSegments() {
		double logHR = 0.0;
		
		List<Integer> removableEdges = new ArrayList<>();
		for (int i = 0; i < networkEdges.size(); i++) {
			NetworkEdge edge = networkEdges.get(i);
			if (!edge.isRootEdge() && edge.childNode.isReassortment()
					&& edge.parentNode.isCoalescence() && edge.hasSegments.cardinality() == 0) {
				removableEdges.add(i);
			}
		}
		
				
		int nrRemoved = 0-removedEdgesOffset.get();

		while (removableEdges.size() > 0) {
			nrRemoved++;

			logHR -= Math.log(1.0 / (removableEdges.size()));
			int edgeInd = Randomizer.nextInt(removableEdges.size());

			logHR += removeEmptyReassortmentEdge(networkEdges.get(removableEdges.get(edgeInd)));

			// In case and invalid network has been created
			if (logHR == Double.NEGATIVE_INFINITY)
				return Double.NEGATIVE_INFINITY;

			removableEdges = new ArrayList<>();
			for (int i = 0; i < networkEdges.size(); i++) {
				NetworkEdge edge = networkEdges.get(i);
				if (!edge.isRootEdge() && edge.childNode.isReassortment()
						&& edge.parentNode.isCoalescence() && edge.hasSegments.cardinality() == 0) {
					removableEdges.add(i);
				}
			}
		}
		
		netChange -= nrRemoved;

		// probability of adding n empty edges in reverse move
		logHR += Math.log(Math.pow(lambda, nrRemoved)) - lambda - logFactorial(nrRemoved);

		return logHR;
	}

	double removeEmptyReassortmentEdge(NetworkEdge edgeToRemove) {
		double logHR = 0.0;

		double sourceTime = edgeToRemove.childNode.getHeight();
		NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
		NetworkEdge destEdge = getSisterEdge(edgeToRemove);
		if (destEdge.childNode == edgeToRemove.childNode)
			destEdge = sourceEdge;
		double destTime = edgeToRemove.parentNode.getHeight();

		// Remove reassortment edge
		logHR += actuallyRemoveEmptyReassortmentEdge(edgeToRemove);

		if (logHR == Double.NEGATIVE_INFINITY)
			return Double.NEGATIVE_INFINITY;

		int nPossibleSourceEdges = networkEdges.size() - 1;

		// change to -1 because after the move, 1 less can be a source edge
		logHR += Math.log(1.0 / (double) nPossibleSourceEdges) + Math.log(1.0 / sourceEdge.getLength());

		int possibleDestEdges = 0;
		for (NetworkEdge e : networkEdges) {
			if (!e.isRootEdge() && e.parentNode.getHeight() >= sourceTime) {
				possibleDestEdges++;
			}
		}
		
		// works
		logHR += Math.log(1.0 / (possibleDestEdges + 1));

		double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

		if (destEdge.isRootEdge()) {
			// works
			logHR += -(1.0 / emptyAlpha) * (destTime - minDestTime) + Math.log(1.0 / emptyAlpha);
		} else {
			logHR += Math.log(1.0 / (destEdge.parentNode.getHeight() - minDestTime));

		}

		return logHR;
	}

	double actuallyRemoveEmptyReassortmentEdge(NetworkEdge edgeToRemove) {

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
		}


//        if (!networkTerminatesAtMRCA()){
//        	
//            return Double.NEGATIVE_INFINITY;
//        }

		return logHR;
	}

	protected boolean networkTerminatesAtMRCA() {
		
		List<NetworkNode> sortedNodes = new ArrayList<>(network.getNodes());
		sortedNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));
		List<NetworkNode> sampleNodes = sortedNodes.stream().filter(NetworkNode::isLeaf).collect(Collectors.toList());
		double maxSampleHeight = sampleNodes.get(sampleNodes.size() - 1).getHeight();

		int lineages = 0;
		for (NetworkNode node : sortedNodes) {
			switch (node.getChildEdges().size()) {
			case 2:
				// Coalescence

				lineages -= 1;
				break;

			case 1:
				// Reassortment

				if (lineages < 2 && node.getHeight() > maxSampleHeight) {
//					System.out.println(node.getHeight());
//					System.out.println(network);
					NetworkEdge newRootEdge = node.getChildEdges().get(0);
					network.setRootEdge(newRootEdge);
//					System.exit(0);
					return false;
				}

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
	
	protected boolean networkTerminatesAtMRCA(NetworkEdge rootEdge) {
		
		if (rootEdge.childNode.getChildEdges().get(0).childNode.isReassortment() &&
				rootEdge.childNode.getChildEdges().get(1).childNode.isReassortment())
			return false;
		return true;

	}


	/**
	 * Check that each edge is ancestral to at least one segment.
	 *
	 * @return true if all edges are ancestral.
	 */
//    public boolean allEdgesAncestral() {
//        Set<NetworkNode> nodeList = networkInput.get().getNodes();
//        for (NetworkNode node : nodeList) {
//            for (NetworkEdge parentEdge : node.getParentEdges()) {
//                if (parentEdge.hasSegments.isEmpty())
//                    return false;
//            }
//        }
//
//        return true;
//    }

//    private int factorial(int k){
//    	int f = 1;
//    	for (int i = 2; i <= k; i++)
//    		f*=k;
//    	return f;
//    }
//    
//    private int alfactorial(int k){
//    	int f = 1;
//    	for (int i = 2; i <= k; i++) {
//    		System.out.println(" .. " +  f + " " + k);
//    		f*=i;
//    	}
//    	return f;
//    }

	private double logFactorial(int k) {
		if (k < 1) {
			return 0.0; // log(0!) = log(1) = 0
		}
		double logF = 0.0;
		for (int i = 1; i <= k; i++) {
			logF += Math.log(i);
		}
		return logF;
	}

}
