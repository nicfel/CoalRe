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
            "Mean of exponential used for choosing root attachment times.",
            Input.Validate.REQUIRED);
    
	public Input<Boolean> localMove = new Input<>("localMove",
			"If true, only reassortment edges that are ancestral to the reassortment event are considered.", false);

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
        	if (localMove.get())
        		logHR = addLocalReassortment();
        	else
        		logHR = addReassortment();
        }else {
        	try {
            	if (localMove.get())
            		logHR = removeLocalReassortment();
            	else
            		logHR = removeReassortment();
			} catch (Exception e) {
				logHR = Double.NEGATIVE_INFINITY;
			}
        }           

        return logHR;
    }
    
    double addReassortment() {
        double logHR = 0.0;


        List<NetworkEdge> possibleSourceEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=2)
                .collect(Collectors.toList());
        

        NetworkEdge sourceEdge = possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size()));
        double sourceTime = Randomizer.nextDouble()*sourceEdge.getLength() + sourceEdge.childNode.getHeight();

        logHR -= Math.log(1.0/(double)possibleSourceEdges.size())
                + Math.log(1.0/sourceEdge.getLength());

        NetworkEdge destEdge = networkEdges.get(Randomizer.nextInt(networkEdges.size()));
        logHR -= Math.log(1.0/networkEdges.size());

        if (!destEdge.isRootEdge() && destEdge.parentNode.getHeight() < sourceTime)
            return Double.NEGATIVE_INFINITY;

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        double destTime;
        if (destEdge.isRootEdge()) {

            destTime = minDestTime + Randomizer.nextExponential(1.0/alpha);
            logHR -= -(1.0/alpha)*(destTime-minDestTime) + Math.log(1.0/alpha);

        } else {

            destTime = Randomizer.nextDouble()*(destEdge.parentNode.getHeight()-minDestTime) + minDestTime;
            logHR -= Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));

        }

        // Create new reassortment edge

        logHR += addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // HR contribution for reverse move
        int nRemovableEdges = (int) networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=1)
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.parentNode.isCoalescence())
                .count();
        logHR += Math.log(1.0/nRemovableEdges);

        return logHR;
    }


    /** same as above, but the reattachment is to ancestral edges only
     * with the goal of minimizing topology changes on the tree
     * @return
     */
    double addLocalReassortment() {
        double logHR = 0.0;


        List<NetworkEdge> possibleSourceEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=2)
                .collect(Collectors.toList());
        

        NetworkEdge sourceEdge = possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size()));
        double sourceTime = Randomizer.nextDouble()*sourceEdge.getLength() + sourceEdge.childNode.getHeight();



        
        logHR -= Math.log(1.0/(double)possibleSourceEdges.size())
                + Math.log(1.0/sourceEdge.getLength());
        
        // make it a loop
    	NetworkEdge destEdge = sourceEdge;

    	double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

    	double destTime;
        destTime = Randomizer.nextDouble()*(destEdge.parentNode.getHeight()-minDestTime) + minDestTime;
        logHR -= Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));


        // Create new reassortment edge

        logHR += addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // HR contribution for reverse move
        int nRemovableEdges = (int) networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.parentNode.isCoalescence())
                .filter(e -> e.hasSegments.cardinality()>=1)
                .filter(e -> getSpouseEdge(e).parentNode == e.parentNode)
                .count();

        logHR += Math.log(1.0/nRemovableEdges);
        
//        System.out.println("removablesdaljd: " + nRemovableEdges + " " + logHR);


        return logHR;
    }

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

    
	double addReassortmentEdge(NetworkEdge sourceEdge, double sourceTime,
                               NetworkEdge destEdge, double destTime) {

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
        BitSet segsToDivert = getRandomConditionedSubset(sourceEdge.hasSegments);
        logHR -= getLogConditionedSubsetProb(sourceEdge.hasSegments);
        logHR += divertSegments(reassortmentEdge, newEdge1, segsToDivert);
        
        networkEdges.add(reassortmentEdge);
        networkEdges.add(newEdge1);
        networkEdges.add(newEdge2);
        		
        		
        		
//        logHR -= addSegmentsToAncestors(reassortmentEdge, segsToDivert);
//        logHR += removeSegmentsFromAncestors(newEdge1, segsToDivert);

        return logHR;
    }

    double removeReassortment() {
        double logHR = 0.0;

        List<NetworkEdge> removableEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=1)
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.parentNode.isCoalescence())
                .collect(Collectors.toList());

        if (removableEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        NetworkEdge edgeToRemove = removableEdges.get(Randomizer.nextInt(removableEdges.size()));
        logHR -= Math.log(1.0/(removableEdges.size()));

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


        int nPossibleSourceEdges = (int) networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=2)
                .count();

        logHR += Math.log(1.0/(double)nPossibleSourceEdges)
                + Math.log(1.0/sourceEdge.getLength());

        logHR += Math.log(1.0/networkEdges.size());

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        if (destEdge.isRootEdge()) {
            logHR += -(1.0/alpha)*(destTime-minDestTime) + Math.log(1.0/alpha);
        } else {
            logHR += Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));
        }

        return logHR;
    }
    
    double removeLocalReassortment() {
        double logHR = 0.0;
        
//        System.out.println(network);

        List<NetworkEdge> removableEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.parentNode.isCoalescence())
                .filter(e -> e.hasSegments.cardinality()>=1)
                .filter(e -> getSpouseEdge(e).parentNode == e.parentNode)
                .collect(Collectors.toList());

        if (removableEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;
        

        NetworkEdge edgeToRemove = removableEdges.get(Randomizer.nextInt(removableEdges.size()));
        logHR -= Math.log(1.0/(removableEdges.size()));

        double sourceTime = edgeToRemove.childNode.getHeight();
        NetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
        NetworkEdge destEdge = sourceEdge;
        
        
//        double destTime = edgeToRemove.parentNode.getHeight();

        // Remove reassortment edge
        logHR += removeReassortmentEdge(edgeToRemove);

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;
        // HR contribution for reverse move


        int nPossibleSourceEdges = (int) networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()>=2)
                .count();

        logHR += Math.log(1.0/(double)nPossibleSourceEdges)
                + Math.log(1.0/sourceEdge.getLength());

//        logHR += Math.log(1.0/networkEdges.size());

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

//        if (destEdge.isRootEdge()) {
//            logHR += -(1.0/alpha)*(destTime-minDestTime) + Math.log(1.0/alpha);
//        } else {
            logHR += Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));
//        }
        
//        System.out.println("removableEdges: " + removableEdges.size() + " " + logHR);


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
//        logHR -= addSegmentsToAncestors(edgeToRemoveSpouse, segsToDivert);
//        logHR += removeSegmentsFromAncestors(edgeToRemove, segsToDivert);
        logHR += getLogConditionedSubsetProb(edgeToRemoveSpouse.hasSegments);

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

