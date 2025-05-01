package coalre.operators;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

public abstract class EmptyEdgesNetworkOperator extends NetworkOperator {

    public Input<Double> emptyAlphaInput = new Input<>("emptyAlpha",
            "Mean of exponential used for choosing root attachment times.",
            0.1);
    
    public Input<Double> lambdaInput = new Input<>("lambda",
            "lambda of the poisson distribution for how many empty edges to add.",
            0.1);
    
    public Input<Boolean> addRemoveEmptyEdgesInput = new Input<>("addRemoveEmptyEdges",
            "adds empty edges before calling the networkproposal and then removes all empty edges at the end again",
            true);

    private double emptyAlpha;
    private double lambda;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        emptyAlpha = emptyAlphaInput.get();
        lambda = lambdaInput.get();
    }

    @Override
    public double proposal() {
    	
		//System.out.println(getID());
       
        double logHR = 0.0;   
        
//        System.out.println("b");
//        System.out.println(network);        


        // Adds empty network edges
        if (addRemoveEmptyEdgesInput.get()){
        	logHR += addEmptyNetworkSegments();
        	
            if (logHR == Double.NEGATIVE_INFINITY)
            	return Double.NEGATIVE_INFINITY;
        }        
         
        // calls the operator
        logHR += networkProposal();
        


        // removes all the empty network edges in the network again
        if (addRemoveEmptyEdgesInput.get()){
            if (logHR == Double.NEGATIVE_INFINITY)
            	return Double.NEGATIVE_INFINITY;
            
        	logHR += RemoveAllEmptyNetworkSegments();
        }
        
        if (logHR == Double.POSITIVE_INFINITY)
        	return Double.NEGATIVE_INFINITY;
        
        // case there are empty edges, which can happen when addRemoveEmptyEdges is false
		if (!allEdgesAncestral()){
            return Double.NEGATIVE_INFINITY;
		}

        if (logHR>Double.NEGATIVE_INFINITY) {
            for (int segIdx=0; segIdx<segmentTrees.size(); segIdx++)
                network.updateSegmentTree(segmentTrees.get(segIdx), segIdx);
        }
        
//        if (logHR>100) {
//        	System.out.println("logHR: " + logHR + " " + this.getID());
//        	System.out.println(network);
//        	System.exit(0);
//        }
        
                		
        return logHR;
    }
    
    private double addEmptyNetworkSegments(){
    	double logHR = 0.0;
    	
    	// randomly sample the number of edges to add
    	int nrEmptyEdges = (int) Randomizer.nextPoisson(lambda);
    	    	
    	for (int i = 0; i < nrEmptyEdges; i ++){
    		logHR += addEmptyReassortment();
    	}  
    	
    	
    	logHR -= Math.log(Math.pow(lambda, nrEmptyEdges)) - lambda - logFactorial(nrEmptyEdges);
    	
    	return logHR;
    }
    

    
    double addEmptyReassortment() {
        double logHR = 0.0;

        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        // add empty reassortment edges to non empty edges
        List<NetworkEdge> possibleSourceEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .collect(Collectors.toList());

        NetworkEdge sourceEdge = possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size()));
        double sourceTime = Randomizer.nextDouble()*sourceEdge.getLength() + sourceEdge.childNode.getHeight();

        logHR -= Math.log(1.0/(double)possibleSourceEdges.size())
                + Math.log(1.0/sourceEdge.getLength());
               
        List<NetworkEdge> possibleDestEdges = networkEdges.stream()
                .collect(Collectors.toList());

        NetworkEdge destEdge = possibleDestEdges.get(Randomizer.nextInt(possibleDestEdges.size()));
    	// works
        logHR -= Math.log(1.0/possibleDestEdges.size());

        
        if (!destEdge.isRootEdge() && destEdge.parentNode.getHeight() < sourceTime)
            return Double.NEGATIVE_INFINITY;

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        double destTime;
        if (destEdge.isRootEdge()) {
        	// works
            destTime = minDestTime + Randomizer.nextExponential(1.0/emptyAlpha);
            logHR -= -(1.0/emptyAlpha)*(destTime-minDestTime) + Math.log(1.0/emptyAlpha);
        } else {
            destTime = Randomizer.nextDouble()*(destEdge.parentNode.getHeight()-minDestTime) + minDestTime;
            logHR -= Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));
        }
        
      

        // Create new reassortment edge
        logHR += addEmptyReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;
        
        // HR contribution for reverse move
        int nRemovableEdges = (int) network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.hasSegments.cardinality()==0)
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.parentNode.isCoalescence())
                .count();

        // works
        logHR += Math.log(1.0/nRemovableEdges);       
        
        if(!networkTerminatesAtMRCA())
        	return Double.NEGATIVE_INFINITY;


        return logHR;
    }
    
    // Only adds the reassortment edge, but does not diverge segments
    double addEmptyReassortmentEdge(NetworkEdge sourceEdge, double sourceTime,
            NetworkEdge destEdge, double destTime) {

		double logHR = 0.0;
		
		network.startEditing(this);
		
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
		
		return logHR;
	}
    
    public double RemoveRandomEmptyNetworkSegments() {
    	double logHR = 0.0;
    	
        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        List<NetworkEdge> removableEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.hasSegments.cardinality()==0)
                .filter(e -> e.parentNode.isCoalescence())
                .collect(Collectors.toList());
        
        if (removableEdges.size()>0){
	        logHR -= Math.log(1.0/(removableEdges.size()));
	        int edgeInd = Randomizer.nextInt(removableEdges.size());            
	    	logHR += removeEmptyReassortmentEdge(removableEdges.get(edgeInd));
	    	networkEdges = new ArrayList<>(network.getEdges());
    	}
	    	
        if (logHR == Double.NEGATIVE_INFINITY)
        	return Double.NEGATIVE_INFINITY;

        
        if(!networkTerminatesAtMRCA())
        	return Double.NEGATIVE_INFINITY;
        
        return logHR;
    }

    
    public double RemoveAllEmptyNetworkSegments() {
    	double logHR = 0.0;
    	
        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        List<NetworkEdge> removableEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.hasSegments.cardinality()==0)
                .filter(e -> e.parentNode.isCoalescence())
                .collect(Collectors.toList());
        
        int nrRemoved = 0;
        
        while (removableEdges.size()>0){
            nrRemoved++;

            logHR -= Math.log(1.0/(removableEdges.size()));
            int edgeInd = Randomizer.nextInt(removableEdges.size());     
            
        	logHR += removeEmptyReassortmentEdge(removableEdges.get(edgeInd));
        	
        	networkEdges = new ArrayList<>(network.getEdges());
        	
        	// In case and invalid network has been created
            if (logHR == Double.NEGATIVE_INFINITY)
            	return Double.NEGATIVE_INFINITY;

            
            removableEdges = networkEdges.stream()
                    .filter(e -> !e.isRootEdge())
                    .filter(e -> e.childNode.isReassortment())
                    .filter(e -> e.hasSegments.cardinality()==0)
                    .filter(e -> e.parentNode.isCoalescence())
                    .collect(Collectors.toList());            
        } 
        
        // probability of adding n empty edges in reverse move
        logHR += Math.log(Math.pow(lambda, nrRemoved)) -lambda -  logFactorial(nrRemoved);
        
//		if (logFactorial(nrRemoved) != Math.log(factorial(nrRemoved))) {
//			System.out.println(nrRemoved + " " + alfactorial(nrRemoved));
//	        System.out.println(logFactorial(nrRemoved) + " " + Math.log(factorial(nrRemoved)));
//		}
//        
//		if (logHR == Double.POSITIVE_INFINITY) {
//			System.out.println(lambda + " " + nrRemoved + " " + ( Math.log(Math.pow(lambda, nrRemoved)) -lambda -  Math.log(factorial(nrRemoved))) );
//			System.out.println(Math.log(Math.pow(lambda, nrRemoved)));
////			System.out.println(altfactorial(nrRemoved));
//			System.out.println(Math.log(factorial(nrRemoved)));
//		}
        
        if (!allEdgesAncestral()){
        	//TODO change to Exception
        	System.err.println("still has empty segments, should not happen ever!");
        	return Double.NEGATIVE_INFINITY;
        }
        
        if(!networkTerminatesAtMRCA())
        	return Double.NEGATIVE_INFINITY;
        
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
        
        // HR contribution for reverse move
        Set<NetworkEdge> finalNetworkEdges = network.getEdges();

        int nPossibleSourceEdges = (int)finalNetworkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .count();
        
        
        // change to -1 because after the move, 1 less can be a source edge
        logHR += Math.log(1.0/(double) nPossibleSourceEdges )
                + Math.log(1.0/sourceEdge.getLength());        

        List<NetworkEdge> possibleDestEdges = finalNetworkEdges.stream()
                .collect(Collectors.toList());

        // works
        logHR += Math.log(1.0/finalNetworkEdges.size());
    

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        if (destEdge.isRootEdge()) {
        	// works
            logHR += -(1.0/emptyAlpha)*(destTime-minDestTime) + Math.log(1.0/emptyAlpha);
        } else {
            logHR += Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));
            
        }

        return logHR;
    }


    double actuallyRemoveEmptyReassortmentEdge(NetworkEdge edgeToRemove) {

        double logHR = 0.0;

        network.startEditing(this);

        NetworkNode nodeToRemove = edgeToRemove.childNode;
        NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
        NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;

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

        if (secondNodeToRemove.getParentEdges().get(0).isRootEdge()) {
            network.setRootEdge(secondEdgeToExtend);

        } else {
            NetworkEdge secondNodeToRemoveParentEdge = secondNodeToRemove.getParentEdges().get(0);
            NetworkNode secondNodeToRemoveParent = secondNodeToRemoveParentEdge.parentNode;
            secondNodeToRemoveParent.removeChildEdge(secondNodeToRemoveParentEdge);
            secondNodeToRemove.removeParentEdge(secondNodeToRemoveParentEdge);

            secondNodeToRemoveParent.addChildEdge(secondEdgeToExtend);
        }

        if (!networkTerminatesAtMRCA()){
            return Double.NEGATIVE_INFINITY;
        }


        return logHR;
    }
 


    protected boolean networkTerminatesAtMRCA() {
        List<NetworkNode> sortedNodes = new ArrayList<>(network.getNodes());
        sortedNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));
        List<NetworkNode> sampleNodes = sortedNodes.stream().filter(NetworkNode::isLeaf).collect(Collectors.toList());
        double maxSampleHeight = sampleNodes.get(sampleNodes.size()-1).getHeight();

        int lineages = 0;
        for (NetworkNode node : sortedNodes) {
            switch(node.getChildEdges().size()) {
                case 2:
                    // Coalescence

                    lineages -= 1;
                    break;

                case 1:
                    // Reassortment

                    if (lineages < 2 && node.getHeight() > maxSampleHeight)
                        return false;

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
    
    /**
     * Check that each edge is ancestral to at least one segment.
     *
     * @return true if all edges are ancestral.
     */
    public boolean allEdgesAncestral() {
        Set<NetworkNode> nodeList = networkInput.get().getNodes();
        for (NetworkNode node : nodeList) {
            for (NetworkEdge parentEdge : node.getParentEdges()) {
                if (parentEdge.hasSegments.isEmpty())
                    return false;
            }
        }

        return true;
    }


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
