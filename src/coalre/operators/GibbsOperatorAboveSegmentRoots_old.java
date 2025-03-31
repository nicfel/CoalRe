package coalre.operators;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.statistics.NetworkStatsLogger;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

public class GibbsOperatorAboveSegmentRoots_old extends NetworkOperator {

    public Input<RealParameter> reassortmentRateInput = new Input<>("reassortmentRate",
            "Rate of reassortment (per lineage per unit time)", Validate.REQUIRED);

    public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
            "Population model to use.", Validate.REQUIRED);
    
    public Input<RealParameter> binomialProbInput = new Input<>("binomialProb",
            "Probability of a given segment choosing a particular parent.");

    public Input<Double> maxHeightRatioInput = new Input<>(
            "maxHeightRatio",
            "if specified, above the ratio, only coalescent events are allowed.", Double.POSITIVE_INFINITY);

    public Input<Double> redFactorInput = new Input<>(
            "redFactor",
            "by how much the recombination rate should be reduced after reaching the maxHeightRatio.", 0.1);



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
        // get the maximum height of the segment tree roots
        double maxHeight = NetworkStatsLogger.getLociMRCA(network);

    	// get all network edges 
        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        // keep only those that coexist at the time of maxHeight
        List<NetworkEdge> startingEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.parentNode.getHeight()>maxHeight)
                .filter(e -> e.childNode.getHeight()<=maxHeight)
               .collect(Collectors.toList());
        
        if (startingEdges.size()==0)
        	return Double.NEGATIVE_INFINITY;
        
       // simulate the rest of the network starting from mxHeight
        double currentTime = maxHeight;
        double timeUntilNextSample = Double.POSITIVE_INFINITY;
        // get the time when the reassortment rates are reduced
        double recChangeTime = maxHeight*maxHeightRatioInput.get();
        double redFactor = 1.0;
        do {
            // get the current propensities
            int k = startingEdges.size();

            double currentTransformedTime = populationFunction.getIntensity(currentTime);
            double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
            double timeToNextCoal = populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;

            
            double timeToNextReass = k>=1 ? Randomizer.nextExponential(k*reassortmentRate.getValue() * redFactor) : Double.POSITIVE_INFINITY;

            // next event time
            double timeUntilNextEvent = Math.min(timeToNextCoal, timeToNextReass);
            if ((timeUntilNextEvent+currentTime)>recChangeTime) {
                currentTime = recChangeTime;
                redFactor *= redFactorInput.get();
                recChangeTime = Double.POSITIVE_INFINITY;
            }else {
                if (timeUntilNextEvent < timeUntilNextSample) {
                    currentTime += timeUntilNextEvent;
                    if (timeUntilNextEvent == timeToNextCoal)
                        coalesce(currentTime, startingEdges);
                    else
                        reassort(currentTime, startingEdges);
                }
            }

        }
        while (startingEdges.size() > 1);
        
        network.setRootEdge(startingEdges.get(0));
//        System.out.println(network.getExtendedNewick());

        
        return Double.POSITIVE_INFINITY;

        

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
        	if (binomialProbInput.get()==null) {
	            if (Randomizer.nextBoolean()) {
	                hasSegs_left.set(segIdx);
	            } else {
	                hasSegs_right.set(segIdx);
	            }
        	}else {
	            if (Randomizer.nextDouble()>binomialProbInput.get().getArrayValue()) {
	                hasSegs_left.set(segIdx);
	            } else {
	                hasSegs_right.set(segIdx);
	            }

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
