package coalre.operators;

import beast.base.core.Input;
import beast.pkgmgmt.Package;
import beast.base.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

public class AddRemoveReassortmentExponential extends DivertSegmentOperator {

	public Input<Double> alphaInput = new Input<>("alpha",
			"Mean of exponential used for choosing root attachment times.", 1.0);

	public Input<Boolean> localMove = new Input<>("localMove",
			"If true, only reassortment edges that are ancestral to the reassortment event are considered.", false);

	public Input<Double> addProbabilityInput = new Input<>("addProbability",
			"Probability of adding a reassortment edge.", 0.5);
	
	private double alpha;

	@Override
	public void initAndValidate() {
		super.initAndValidate();

		alpha = alphaInput.get();
	}

	@Override
	public double networkProposal() {

		double logHR;
//		System.out.println(network);
		if (Randomizer.nextDouble() < addProbabilityInput.get()) {
				logHR = addReassortment();
		} else {
				logHR = removeReassortment();
//			} catch (Exception e) {
//				// dangerous implementation!!
//				logHR = Double.NEGATIVE_INFINITY;
//			}
		}
//		System.out.println(network);

		return logHR;
	}

	double addReassortment() {
		double logHR = 0.0;

		List<NetworkEdge> possibleSourceEdges = networkEdges.stream().filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() >= 2).collect(Collectors.toList());

		NetworkEdge sourceEdge = possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size()));
		
		logHR -= Math.log(1.0 / (double) possibleSourceEdges.size());

		
		// sample two times t1 and t2
		double t1 =  Randomizer.nextDouble() * sourceEdge.getLength();
		double t2 =  Randomizer.nextDouble() * sourceEdge.getLength();
		
		logHR -= Math.log(2.0 / (sourceEdge.getLength() * sourceEdge.getLength()));
		
		double sourceTime = Math.min(t1,t2) + sourceEdge.childNode.getHeight();
		double destTime = Math.max(t1,t2) + sourceEdge.childNode.getHeight();	
		
		

		List<NetworkEdge> destEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.getHeight() > destTime)
				.filter(e -> e.childNode.getHeight() < destTime)
				.collect(Collectors.toList());	

		
		NetworkEdge destEdge = destEdges.get(Randomizer.nextInt(destEdges.size()));
		logHR -= Math.log(1.0 / destEdges.size());

		// Create new reassortment edge
		logHR += addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);

		if (logHR == Double.NEGATIVE_INFINITY)
			return Double.NEGATIVE_INFINITY;

		// HR contribution for reverse move
		int nRemovableEdges;
		
		nRemovableEdges = (int) networkEdges.stream().filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() >= 1)
				.filter(e -> e.childNode.isReassortment())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> getSpouseEdge(e).parentNode.getHeight() >= e.parentNode.getHeight()) //= to allow for loops
				.count();

		logHR += Math.log(1.0 / nRemovableEdges);

		return logHR;
	}

	/**
	 * same as above, but the reattachment is to ancestral edges only with the goal
	 * of minimizing topology changes on the tree
	 * 
	 * @return
	 */
	private void getAllAncestralEdges(NetworkEdge edge, List<NetworkEdge> destEdges) {
		if (destEdges.contains(edge)) {
			return;
		}
		destEdges.add(edge);
		if (edge.parentNode != null) {
			for (NetworkEdge parentEdge : edge.parentNode.getParentEdges()) {
				getAllAncestralEdges(parentEdge, destEdges);
			}
		}

	}

	double addReassortmentEdge(NetworkEdge sourceEdge, double sourceTime, NetworkEdge destEdge, double destTime) {

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

        // Choose segments to divert to new edge
        if (divertOneSegmentInput.get()) {
        	double prob = 1./sourceEdge.hasSegments.cardinality();
            BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments, prob);
            logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments, segsToDivert, prob);
            logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
        }else {
	        BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
            logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments);
            logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
        }

		networkEdges.add(reassortmentEdge);
		networkEdges.add(newEdge1);
		networkEdges.add(newEdge2);

		return logHR;
	}

	double removeReassortment() {
		double logHR = 0.0;

		
        List<NetworkEdge> removableEdges = networkEdges.stream()
        		.filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() >= 1)
				.filter(e -> e.childNode.isReassortment())
				.filter(e -> e.parentNode.isCoalescence())
				.filter(e -> getSpouseEdge(e).parentNode.getHeight() >= e.parentNode.getHeight())
				.collect(Collectors.toList());

        
		

		if (removableEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;

		NetworkEdge edgeToRemove = removableEdges.get(Randomizer.nextInt(removableEdges.size()));

		logHR -= Math.log(1.0 / (removableEdges.size()));

		double destTime = edgeToRemove.parentNode.getHeight();
		NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
		NetworkEdge destEdge = getSisterEdge(edgeToRemove);
		if (destEdge.childNode == edgeToRemove.childNode)
			destEdge = sourceEdge;

		// Remove reassortment edge
		logHR += removeReassortmentEdge(edgeToRemove);

		if (logHR == Double.NEGATIVE_INFINITY)
			return Double.NEGATIVE_INFINITY;

		// HR contribution for reverse move

		int nPossibleSourceEdges = (int) networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() >= 2).count();

		logHR += Math.log(1.0 / (double) nPossibleSourceEdges);

		
		logHR += Math.log(2.0 / (sourceEdge.getLength() * sourceEdge.getLength()));

		int destEdgesNumber = (int) networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.getHeight() > destTime)
				.filter(e -> e.childNode.getHeight() < destTime)
				.count();	
		
		
		logHR += Math.log(1.0 / destEdgesNumber);


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

//        if (!networkTerminatesAtMRCA())
//            return Double.NEGATIVE_INFINITY;

		return logHR;
	}
	
//	double combineReassortment() {
//		double logHR = 0.0;
//		
//		List<NetworkEdge> potentialSourceEdges = networkEdges.stream()
//				.filter(e -> !e.isRootEdge())
//				.filter(e -> e.parentNode.isReassortment())
//				.filter(e -> e.childNode.isReassortment())
//				.collect(Collectors.toList());
//		
//		if (potentialSourceEdges.isEmpty())
//			return Double.NEGATIVE_INFINITY;
//		
//		
//		
//		
//		
//		
//		NetworkEdge sourceEdge = potentialSourceEdges.get(Randomizer.nextInt(potentialSourceEdges.size()));
//		logHR -= Math.log(1.0 / potentialSourceEdges.size());
//					
//		NetworkEdge destEdge = sourceEdge;
//		
//		double minTime = Math.max(newSourceTime, destEdge.childNode.getHeight());
//		
//		double newDestTime = Randomizer.nextDouble() * (destEdge.parentNode.getHeight()-minTime) + minTime;
//		logHR -= Math.log(1.0 / (destEdge.parentNode.getHeight()-minTime));		
//		
//
//		// Create new reassortment edge
//		logHR += addReassortmentEdge(sourceEdge, newSourceTime, destEdge, newDestTime);
//		if (logHR == Double.NEGATIVE_INFINITY)
//			return Double.NEGATIVE_INFINITY;
//		
//		int nRemovableEdges = (int) networkEdges.stream()
//				.filter(e -> !e.isRootEdge())
//				.filter(e -> e.childNode.isReassortment())
//				.filter(e -> getSpouseEdge(e)==getSisterEdge(e))
//				.filter(e -> e.hasSegments.cardinality() >= 1)
//				.filter(e -> getSpouseEdge(e).hasSegments.cardinality() >= 1)
//				.count();
//		// write the same into a for loop
//
//		logHR += Math.log(1.0 / nRemovableEdges);
//
//		return logHR;
//	}
//
//	double splitReassortment() {
//		double logHR = 0.0;
//
//		List<NetworkEdge> removableEdges = networkEdges.stream()
//				.filter(e -> !e.isRootEdge())
//				.filter(e -> e.parentNode.isCoalescence())				
//				.filter(e -> e.childNode.isReassortment())
//				.filter(e -> getSpouseEdge(e)==getSisterEdge(e))
//				.filter(e -> e.hasSegments.cardinality() >= 1)
//				.filter(e -> getSpouseEdge(e).hasSegments.cardinality() >= 1)
//				.collect(Collectors.toList());
//
//		if (removableEdges.isEmpty())
//			return Double.NEGATIVE_INFINITY;
//		
//		logHR -= Math.log(1.0 / (removableEdges.size()));
//		
//
//		NetworkEdge edgeToRemove = removableEdges.get(Randomizer.nextInt(removableEdges.size()));		
//		NetworkEdge childEdge = edgeToRemove.childNode.getChildEdges().get(0);
//		NetworkEdge newDestEdge = getSisterEdge(edgeToRemove);
//		
//		double oldSourceTime = edgeToRemove.childNode.getHeight();
//		double oldParentHeight = edgeToRemove.parentNode.getParentEdges().get(0).parentNode.getHeight();
//				
//		logHR += removeReassortmentEdge(edgeToRemove);
//
//		if (logHR == Double.NEGATIVE_INFINITY)
//			return Double.NEGATIVE_INFINITY;
//		
//
//
//		// HR contribution for reverse move
//		logHR += Math.log(1.0 / (oldParentHeight - oldSourceTime));	
//		
//		logHR += Math.log(1.0 / (network.getRootEdge().childNode.getHeight()));
//
//		
//		// get the possible desEdges for reverse calculation
//		int nPossibleDestEdges = (int) networkEdges.stream()
//				.filter(e -> !e.isRootEdge())
//				.filter(e -> e.childNode.getHeight() < oldSourceTime)
//				.filter(e -> e.parentNode.getHeight() > oldSourceTime)
//				.filter(e -> e.hasSegments.cardinality() > 1)
//                .count();
//		
//		logHR += Math.log(1.0 / nPossibleDestEdges);
//
//		return logHR;
//	}
//
//	
//	
	
//	
//	double addLocalReassortmentEdge(NetworkEdge sourceEdge, double sourceTime, NetworkEdge destEdge, double destTime, int segment) {
//
//		double logHR = 0.0;
//
//		NetworkNode sourceNode = new NetworkNode();
//		sourceNode.setHeight(sourceTime);
//
//		NetworkNode oldSourceEdgeParent = sourceEdge.parentNode;
//		oldSourceEdgeParent.removeChildEdge(sourceEdge);
//		sourceNode.addChildEdge(sourceEdge);
//
//		NetworkEdge newEdge1 = new NetworkEdge();
//		sourceNode.addParentEdge(newEdge1);
//		oldSourceEdgeParent.addChildEdge(newEdge1);
//
//		newEdge1.hasSegments = (BitSet) sourceEdge.hasSegments.clone();
//
//		if (destEdge == sourceEdge)
//			destEdge = newEdge1;
//
//		NetworkNode destNode = new NetworkNode();
//		destNode.setHeight(destTime);
//
//		NetworkNode oldDestEdgeParent = destEdge.parentNode;
//		if (oldDestEdgeParent != null) {
//			oldDestEdgeParent.removeChildEdge(destEdge);
//		}
//
//		destNode.addChildEdge(destEdge);
//
//		NetworkEdge newEdge2 = new NetworkEdge();
//		destNode.addParentEdge(newEdge2);
//
//		if (oldDestEdgeParent == null) {
//			network.setRootEdge(newEdge2);
//		} else {
//			oldDestEdgeParent.addChildEdge(newEdge2);
//		}
//
//		newEdge2.hasSegments = (BitSet) destEdge.hasSegments.clone();
//
//		NetworkEdge reassortmentEdge = new NetworkEdge();
//		sourceNode.addParentEdge(reassortmentEdge);
//		destNode.addChildEdge(reassortmentEdge);
//		reassortmentEdge.hasSegments = new BitSet();
//		
////		BitSet segsToDivert = new BitSet();
//
//		// ensure this segment gets diverted
//		BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
//		logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments);
//		
//		
//		
//		
//		logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
//
//		
//		
////		segsToDivert.set(segment);
////		logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
//		
//
//		networkEdges.add(reassortmentEdge);
//		networkEdges.add(newEdge1);
//		networkEdges.add(newEdge2);
//
//		return logHR;
//	}
//	
//	double removeLocalReassortmentEdge(NetworkEdge edgeToRemove) {
//		double logHR = 0.0;
//
//		NetworkNode nodeToRemove = edgeToRemove.childNode;
//		NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
//		NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;
//		networkEdges.remove(edgeToRemove);
//		networkEdges.remove(edgeToRemoveSpouse);
//
//		// Divert segments away from chosen edge
//		BitSet segsToDivert = (BitSet) edgeToRemove.hasSegments.clone();
//		logHR += divertSegments(edgeToRemoveSpouse, edgeToRemove, segsToDivert);
//		logHR += getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments);
//
//		// Remove edge and associated nodes
//		NetworkEdge edgeToExtend = nodeToRemove.getChildEdges().get(0);
//		nodeToRemove.removeChildEdge(edgeToExtend);
//		nodeToRemove.removeParentEdge(edgeToRemove);
//		nodeToRemove.removeParentEdge(edgeToRemoveSpouse);
//		edgeToRemoveSpouseParent.removeChildEdge(edgeToRemoveSpouse);
//		edgeToRemoveSpouseParent.addChildEdge(edgeToExtend);
//
//		NetworkNode secondNodeToRemove = edgeToRemove.parentNode;
//		NetworkEdge secondEdgeToExtend = getSisterEdge(edgeToRemove);
//
//		secondNodeToRemove.removeChildEdge(secondEdgeToExtend);
//		secondNodeToRemove.removeChildEdge(edgeToRemove);
//
//		networkEdges.remove(secondNodeToRemove.getParentEdges().get(0));
//		if (secondNodeToRemove.getParentEdges().get(0).isRootEdge()) {
//			network.setRootEdge(secondEdgeToExtend);
//
//		} else {
//			NetworkEdge secondNodeToRemoveParentEdge = secondNodeToRemove.getParentEdges().get(0);
//			NetworkNode secondNodeToRemoveParent = secondNodeToRemoveParentEdge.parentNode;
//			secondNodeToRemoveParent.removeChildEdge(secondNodeToRemoveParentEdge);
//			secondNodeToRemove.removeParentEdge(secondNodeToRemoveParentEdge);
//
//			secondNodeToRemoveParent.addChildEdge(secondEdgeToExtend);
//			networkEdges.remove(secondNodeToRemoveParentEdge);
//		}
//
////        if (!networkTerminatesAtMRCA())
////            return Double.NEGATIVE_INFINITY;
//
//		return logHR;
//	}
////	
//	boolean validLocalReassortmentEdge(NetworkEdge edge) {
//		NetworkEdge spouseEdge = getSpouseEdge(edge);
//		
//		List<NetworkEdge> destEdges = new ArrayList<>();
//		if (spouseEdge.parentNode.isCoalescence())
//			destEdges.add(getSisterEdge(edge));
//		
//		for (NetworkEdge parentEdge : spouseEdge.parentNode.getParentEdges()) {
//			if (parentEdge.isRootEdge())
//				continue;
//			destEdges.add(parentEdge);
//
//			if (parentEdge.parentNode.isCoalescence())
//				destEdges.add(getSisterEdge(parentEdge));
//		}
//		
//		if (destEdges.isEmpty())
//			return false;
//		
//			
//		return destEdges.contains(edge.parentNode.getParentEdges().get(0));
//	}
}
