package coalre.operators;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Package;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

public class GibbsOperatorAboveSegmentRoots extends NetworkOperator {

    public Input<RealParameter> reassortmentRateInput = new Input<>("reassortmentRate",
            "Rate of reassortment (per lineage per unit time)", Validate.REQUIRED);

    public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
            "Population model to use.", Validate.REQUIRED);

    private int nSegments;
    
    private PopulationFunction populationFunction;
    private RealParameter reassortmentRate;

    @Override
    public void initAndValidate() {
    	nSegments = segmentTreesInput.get().size();
    	
        populationFunction = populationFunctionInput.get();
        reassortmentRate = reassortmentRateInput.get();

    	
        super.initAndValidate();
    }

    @Override
    public double networkProposal() {
    	return resimulate();
    	
    }

    double resimulate() {
    	network.startEditing(this);
    	
    	// get the place where to cut
    	double maxHeight = getMaxSegmentMRCA();

    	// get all network edges 
        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        // keep only those that coexist at the time of maxHeight
        List<NetworkEdge> startingEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.getHeight()>maxHeight)
                .collect(Collectors.toList());
        
        if (startingEdges.size()==0)
        	return Double.NEGATIVE_INFINITY;
                
        // simulate the rest of the network starting from mxHeight
        double currentTime = maxHeight;
        double timeUntilNextSample = Double.POSITIVE_INFINITY;
        do {

            // get the current propensities
            int k = startingEdges.size();

            double currentTransformedTime = populationFunction.getIntensity(currentTime);
            double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
            double timeToNextCoal = populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;

            double timeToNextReass = k>=1 ? Randomizer.nextExponential(k*reassortmentRate.getValue()) : Double.POSITIVE_INFINITY;

            // next event time
            double timeUntilNextEvent = Math.min(timeToNextCoal, timeToNextReass);
            if (timeUntilNextEvent < timeUntilNextSample) {
                currentTime += timeUntilNextEvent;
                if (timeUntilNextEvent == timeToNextCoal)
                    coalesce(currentTime, startingEdges);
                else
                    reassort(currentTime, startingEdges);
            }

        }
        while (startingEdges.size() > 1);
        
        network.setRootEdge(startingEdges.get(0));

        
        return Double.POSITIVE_INFINITY;

        

    }

    double getMaxSegmentMRCA(){
    	double maxHeight = 0.0;
    	for (int i = 0; i < segmentTreesInput.get().size(); i++){
    		double height = segmentTreesInput.get().get(i).getRoot().getHeight();
    		if (height>maxHeight)
    			maxHeight=height;
    	}
    	
    	return maxHeight;
    }
    
    private void coalesce(double coalescentTime, List<NetworkEdge> extantLineages) {
        // Sample the pair of lineages that are coalescing:
        NetworkEdge lineage1 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
        NetworkEdge lineage2;
        do {
            lineage2 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
        } while (lineage1 == lineage2);

        // Create coalescent node
        NetworkNode coalescentNode = new NetworkNode();
        coalescentNode.setHeight(coalescentTime)
                .addChildEdge(lineage1)
                .addChildEdge(lineage2);
        lineage1.parentNode = coalescentNode;
        lineage2.parentNode = coalescentNode;

        // Merge segment flags:
        BitSet hasSegments = new BitSet();
        hasSegments.or(lineage1.hasSegments);
        hasSegments.or(lineage2.hasSegments);

        // Create new lineage
        NetworkEdge lineage = new NetworkEdge(null, coalescentNode, hasSegments);
        coalescentNode.addParentEdge(lineage);

        extantLineages.remove(lineage1);
        extantLineages.remove(lineage2);
        extantLineages.add(lineage);
    }
    
    private void sample(List<NetworkNode> remainingSampleNodes, List<NetworkEdge> extantLineages) {
        // sample the network node
        NetworkNode n = remainingSampleNodes.get(0);

        // Create corresponding lineage
        BitSet hasSegs = new BitSet();
        hasSegs.set(0, nSegments);
        NetworkEdge lineage = new NetworkEdge(null, n, hasSegs);
        extantLineages.add(lineage);
        n.addParentEdge(lineage);

        remainingSampleNodes.remove(0);
    }


    private void reassort(double reassortmentTime, List<NetworkEdge> extantLineages) {
        NetworkEdge lineage = extantLineages.get(Randomizer.nextInt(extantLineages.size()));

        BitSet hasSegs_left = new BitSet();
        BitSet hasSegs_right = new BitSet();

        for (int segIdx = lineage.hasSegments.nextSetBit(0);
             segIdx != -1; segIdx = lineage.hasSegments.nextSetBit(segIdx+1)) {
            if (Randomizer.nextBoolean()) {
                hasSegs_left.set(segIdx);
            } else {
                hasSegs_right.set(segIdx);
            }
        }

        // Stop here if reassortment event is unobservable
        if (hasSegs_left.cardinality() == 0 || hasSegs_right.cardinality() == 0)
            return;

        // Create reassortment node
        NetworkNode node = new NetworkNode();
        node.setHeight(reassortmentTime).addChildEdge(lineage);

        // Create reassortment lineages
        NetworkEdge leftLineage = new NetworkEdge(null, node, hasSegs_left);
        NetworkEdge rightLineage = new NetworkEdge(null, node, hasSegs_right);
        node.addParentEdge(leftLineage);
        node.addParentEdge(rightLineage);

        extantLineages.remove(lineage);
        extantLineages.add(leftLineage);
        extantLineages.add(rightLineage);
    }
    
}
