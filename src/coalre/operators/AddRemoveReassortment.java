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

public class AddRemoveReassortment extends DivertSegmentOperator {

	public Input<Double> alphaInput = new Input<>("alpha",
			"Mean of exponential used for choosing root attachment times.", 1.0);

	public Input<Boolean> localMove = new Input<>("localMove",
			"If true, only reassortment edges that are ancestral to the reassortment event are considered.", false);

	public Input<Double> addProbabilityInput = new Input<>("addProbability",
			"Probability of adding a reassortment edge.", 0.5);
	
	public Input<Integer> edgeDistanceInput = new Input<>("edgeDistance",
			"Number of edges to traverse when selecting destination edge for local move.", 1);
	
	private double alpha;

	@Override
	public void initAndValidate() {
		super.initAndValidate();

		alpha = alphaInput.get();
	}

	@Override
	public double networkProposal() {

		double logHR;
		if (Randomizer.nextBoolean()) {
			if (localMove.get()) {
//				System.out.println(network);
				logHR = addLocalReassortment();
//				System.out.println(network+"\n");
			}else {
				logHR = addReassortment();
			}
		} else {
			if (localMove.get())
				logHR = removeLocalReassortment();
			else
				logHR = removeReassortment();
		}

		return logHR;
	}

	double addReassortment() {
		double logHR = 0.0;

		List<NetworkEdge> possibleSourceEdges = networkEdges.stream().filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() >= 2).collect(Collectors.toList());

		NetworkEdge sourceEdge = possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size()));
		double sourceTime = Randomizer.nextDouble() * sourceEdge.getLength() + sourceEdge.childNode.getHeight();

		logHR -= Math.log(1.0 / (double) possibleSourceEdges.size()) + Math.log(1.0 / sourceEdge.getLength());

		List<NetworkEdge> destEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.getHeight() >= sourceTime)
				.collect(Collectors.toList());	
		destEdges.add(network.getRootEdge());
		
		NetworkEdge destEdge = destEdges.get(Randomizer.nextInt(destEdges.size()));
		logHR -= Math.log(1.0 / destEdges.size());


		double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

		double destTime;
		if (destEdge.isRootEdge()) {
			destTime = minDestTime + Randomizer.nextExponential(1.0 / alpha);
			logHR -= -(1.0 / alpha) * (destTime - minDestTime) + Math.log(1.0 / alpha);
		} else {
			destTime = Randomizer.nextDouble() * (destEdge.parentNode.getHeight() - minDestTime) + minDestTime;
			logHR -= Math.log(1.0 / (destEdge.parentNode.getHeight() - minDestTime));
		}

		// Create new reassortment edge
		logHR += addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);

		if (logHR == Double.NEGATIVE_INFINITY)
			return Double.NEGATIVE_INFINITY;

		// HR contribution for reverse move
		int nRemovableEdges;
		
		nRemovableEdges= (int) networkEdges.stream().filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() >= 1).filter(e -> e.childNode.isReassortment())
				.filter(e -> e.parentNode.isCoalescence()).count();

		logHR += Math.log(1.0 / nRemovableEdges);

		return logHR;
	}


	double removeReassortment() {
		double logHR = 0.0;

		
		List<NetworkEdge> removableEdges;
        removableEdges = networkEdges.stream().filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() >= 1).filter(e -> e.childNode.isReassortment())
				.filter(e -> e.parentNode.isCoalescence()).collect(Collectors.toList());


		if (removableEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;

		NetworkEdge edgeToRemove = removableEdges.get(Randomizer.nextInt(removableEdges.size()));

		logHR -= Math.log(1.0 / (removableEdges.size()));

		double sourceTime = edgeToRemove.childNode.getHeight();
		NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
		NetworkEdge destEdge = getSisterEdge(edgeToRemove);
		if (destEdge.childNode == edgeToRemove.childNode)
			destEdge = sourceEdge;
		double destTime = edgeToRemove.parentNode.getHeight();

		// Remove reassortment edge
		logHR += removeReassortmentEdge(edgeToRemove);

		if (logHR == Double.NEGATIVE_INFINITY)
			return Double.NEGATIVE_INFINITY;

		// HR contribution for reverse move

		int nPossibleSourceEdges = (int) networkEdges.stream().filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() >= 2).count();

		logHR += Math.log(1.0 / (double) nPossibleSourceEdges) + Math.log(1.0 / sourceEdge.getLength());

		

		int destEdgesNumber = (int) networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.parentNode.getHeight() >= sourceTime)
				.count();	
		
		destEdgesNumber += 1; // include root edge
		
		logHR += Math.log(1.0 / destEdgesNumber);

		double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

		if (destEdge.isRootEdge()) {
			logHR += -(1.0 / alpha) * (destTime - minDestTime) + Math.log(1.0 / alpha);
		} else {
			logHR += Math.log(1.0 / (destEdge.parentNode.getHeight() - minDestTime));
		}

		return logHR;
	}

	double addLocalReassortment() {
		double logHR = 0.0;
		
		
		// 1. get all possible source edges		
		List<NetworkEdge> potentialSourceEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() > 1)		
				.collect(Collectors.toList());
		
		// if no source edges, return negative infinity
		if (potentialSourceEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;
		
		// select source edge and count HR contribution
		NetworkEdge sourceEdge = potentialSourceEdges.get(Randomizer.nextInt(potentialSourceEdges.size()));
		logHR -= Math.log(1.0 / potentialSourceEdges.size());
		
		// 2.select source time and count HR contribution
		double newSourceTime = Randomizer.nextDouble() * sourceEdge.getLength() + sourceEdge.childNode.getHeight();
		logHR -= Math.log(1.0 / sourceEdge.getLength());

		// get all target edges
		List<NetworkEdge> potentialDestEdges = new ArrayList<>();
//		for (NetworkEdge e : sourceEdge.parentNode.getParentEdges())
		getTargetEdgesUp(sourceEdge, potentialDestEdges, edgeDistanceInput.get(), newSourceTime);
        // only keep unique target Edges
		potentialDestEdges = potentialDestEdges.stream()
                .distinct()
                .collect(Collectors.toList());

		if (potentialDestEdges.isEmpty())
			return Double.NEGATIVE_INFINITY;

		// 3. select dest edge
        NetworkEdge destEdge = potentialDestEdges.get(Randomizer.nextInt(potentialDestEdges.size()));
        logHR -= Math.log(1.0/potentialDestEdges.size());

        // 4. select dest time
        double minTime = Math.max(newSourceTime, destEdge.childNode.getHeight());
		double newDestTime;
		if (destEdge.isRootEdge()) {
			newDestTime = minTime + Randomizer.nextExponential(1.0 / alpha);
			logHR -= -(1.0 / alpha) * (newDestTime - minTime) + Math.log(1.0 / alpha);
			
		} else {
			newDestTime = Randomizer.nextDouble() * (destEdge.parentNode.getHeight()-minTime) + minTime;
			logHR -= Math.log(1.0 / (destEdge.parentNode.getHeight()-minTime));
		}		
		// 5. Create new reassortment edge
		
		logHR += addReassortmentEdge(sourceEdge, newSourceTime, destEdge, newDestTime);
		if (logHR == Double.NEGATIVE_INFINITY)
			return Double.NEGATIVE_INFINITY;
				
		// 6. Select removable edge - only those that could have been added by deltaAddRecombination
        List<Integer> removableEdges = new ArrayList<>();
        for (int i = 0; i < networkEdges.size(); i++) {
            NetworkEdge edge = networkEdges.get(i);
            if (!edge.isRootEdge() && edge.childNode.isReassortment()
                    && edge.parentNode.isCoalescence() && edge.hasSegments.cardinality() >= 1) {

                NetworkEdge potentialSourceEdge = getSpouseEdge(edge);
                NetworkEdge trueDestEdge = getSisterEdge(edge);

                // Check if destEdge is reachable from sourceEdge via delta traversal
                List<NetworkEdge> reachableEdges = new ArrayList<>();
                int dist = edgeDistanceInput.get();
        			getTargetEdgesUp(potentialSourceEdge, reachableEdges, dist, edge.childNode.getHeight());

                boolean isReachable = reachableEdges.contains(trueDestEdge);
                if (!isReachable && !trueDestEdge.isRootEdge()) {
                    isReachable = reachableEdges.contains(trueDestEdge.parentNode.getParentEdges().get(0));
                }
                if (isReachable) {
                    removableEdges.add(i);
                }
            }
        }
        
        
        // 7. HR contribution for reverse move for the amount of removable edges
        logHR += Math.log(1.0/removableEdges.size());

        return logHR;
	}

	double removeLocalReassortment() {
		double logHR = 0.0;

		// 1. get all possible removable edges
        List<Integer> removableEdges = new ArrayList<>();
        for (int i = 0; i < networkEdges.size(); i++) {
            NetworkEdge edge = networkEdges.get(i);
            if (!edge.isRootEdge() && edge.childNode.isReassortment()
                    && edge.parentNode.isCoalescence() && edge.hasSegments.cardinality() >= 1) {

                NetworkEdge potentialSourceEdge = getSpouseEdge(edge);
                NetworkEdge trueDestEdge = getSisterEdge(edge);

                // 2. Check if destEdge is reachable from sourceEdge via delta traversal
                List<NetworkEdge> reachableEdges = new ArrayList<>();
                int dist = edgeDistanceInput.get();
        			getTargetEdgesUp(potentialSourceEdge, reachableEdges, dist, edge.childNode.getHeight());

                boolean isReachable = reachableEdges.contains(trueDestEdge);
                if (!isReachable && !trueDestEdge.isRootEdge()) {
                    isReachable = reachableEdges.contains(trueDestEdge.parentNode.getParentEdges().get(0));
                }
                if (isReachable) {
                    removableEdges.add(i);
                }
            }
        }


        // only keep unique edges
        if (removableEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        NetworkEdge edgeToRemove = networkEdges.get(removableEdges.get(Randomizer.nextInt(removableEdges.size())));
        logHR -= Math.log(1.0/(removableEdges.size()));



		// 3. Extract source/dest information
        double sourceTime = edgeToRemove.childNode.getHeight();
        NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
        NetworkEdge destEdge = getSisterEdge(edgeToRemove);
        if (destEdge.childNode == edgeToRemove.childNode)
            destEdge = sourceEdge;
        double destTime = edgeToRemove.parentNode.getHeight();

		logHR += removeReassortmentEdge(edgeToRemove);

		if (logHR == Double.NEGATIVE_INFINITY)
			return Double.NEGATIVE_INFINITY;


		// 4. get all possible reverse target edges
        List<NetworkEdge> reverseTargetEdges = new ArrayList<>();
    	getTargetEdgesUp(sourceEdge, reverseTargetEdges, edgeDistanceInput.get(), sourceTime);
        reverseTargetEdges = reverseTargetEdges.stream()
                .distinct()
                .collect(Collectors.toList());
		// 6. HR contribution for reverse move for the amount of reverse target edges
        logHR += Math.log(1.0 / reverseTargetEdges.size());
		// 7. HR contribution for reverse move for the amount of source edge length
        logHR += Math.log(1.0 / sourceEdge.getLength());
		// 8. HR contribution for reverse move for the amount of dest edge height
        double minTime = Math.max(destEdge.childNode.getHeight(), sourceTime);
		if (destEdge.isRootEdge()) {
			logHR += -(1.0 / alpha) * (destTime - minTime) + Math.log(1.0 / alpha);
		} else {
			logHR += Math.log(1.0 / (destEdge.parentNode.getHeight()-minTime));
		}

		// 9. get the possible desEdges for reverse calculation
		int nPossibleDestEdges = (int) networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> e.hasSegments.cardinality() > 1)
                .count();
		logHR += Math.log(1.0 / nPossibleDestEdges);

		return logHR;
	}

	
	
    // Traverse up the network (forward in time) - only direct ancestors
    private void getTargetEdgesUp(NetworkEdge currentEdge,
                                List<NetworkEdge> targetEdges, int remainingEdgeCount, double minHeight) {

        remainingEdgeCount--;

        // Check if we've hit a root edge (parentNode is null)
        if (currentEdge.isRootEdge()) {
            // Include root edge as potential target
            targetEdges.add(currentEdge);
            return;
        }

        if (remainingEdgeCount < 0 || currentEdge.parentNode.getHeight() < minHeight) {
            return;
        }

        targetEdges.add(currentEdge);

		if (currentEdge.parentNode.isCoalescence()) {
			NetworkEdge sisterEdge = getSisterEdge(currentEdge);
			getTargetEdgesDown(sisterEdge, targetEdges, remainingEdgeCount, minHeight);
			getTargetEdgesUp(currentEdge.parentNode.getParentEdges().get(0), targetEdges, remainingEdgeCount, minHeight);
		}else {
			getTargetEdgesUp(currentEdge.parentNode.getParentEdges().get(0), targetEdges, remainingEdgeCount, minHeight);
			getTargetEdgesUp(currentEdge.parentNode.getParentEdges().get(1), targetEdges, remainingEdgeCount, minHeight);
		}

    }
    
	private void getTargetEdgesDown(NetworkEdge currentEdge,
            List<NetworkEdge> targetEdges, int remainingEdgeCount, double minHeight) {
        remainingEdgeCount--;

        if (currentEdge.isRootEdge() || remainingEdgeCount < 0 || currentEdge.parentNode.getHeight() < minHeight) {
            return;
        }

        targetEdges.add(currentEdge);
        if (currentEdge.childNode.isCoalescence()) {
			getTargetEdgesDown(currentEdge.childNode.getChildEdges().get(0), targetEdges, remainingEdgeCount, minHeight);
			getTargetEdgesDown(currentEdge.childNode.getChildEdges().get(1), targetEdges, remainingEdgeCount, minHeight);
        }else if (currentEdge.childNode.isReassortment()) {
        	getTargetEdgesDown(currentEdge.childNode.getChildEdges().get(0), targetEdges, remainingEdgeCount, minHeight);
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
			BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments, 1.0/sourceEdge.hasSegments.cardinality());
			logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments, segsToDivert, 1.0/ sourceEdge.hasSegments.cardinality());
			logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
		} else {
			BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
			logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments);
			logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
		}

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
		logHR += divertSegments(edgeToRemoveSpouse, edgeToRemove, segsToDivert);
		if (divertOneSegmentInput.get()) {
//			logHR += Math.log(1.0 / edgeToRemoveSpouse.hasSegments.cardinality());
			logHR += getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments, segsToDivert, 1.0/edgeToRemoveSpouse.hasSegments.cardinality());
		} else {
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
	
	
}
