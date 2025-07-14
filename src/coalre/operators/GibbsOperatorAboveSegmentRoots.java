package coalre.operators;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import coalre.distribution.CoalescentWithReassortment;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.statistics.NetworkStatsLogger;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

public class GibbsOperatorAboveSegmentRoots extends NetworkOperator {

    public Input<CoalescentWithReassortment> coalescentDistrInput = new Input<>("coalescentWithReassortment",
            "Mean of exponential used for choosing root attachment times.",
            Input.Validate.REQUIRED);
    

    
    private PopulationFunction populationFunction;
    private Function reassortmentRate;
    public PopulationFunction timeVaryingReassortmentRates;
    
    private boolean isTimeVarying = false;
    private Function binomialProb;
    private double redFactor;
    private double maxHeightRatio;
    private double maxHeight;


    @Override
    public void initAndValidate() {
    	
        populationFunction = coalescentDistrInput.get().populationFunctionInput.get();
        if (coalescentDistrInput.get().timeVaryingReassortmentRatesInput.get()!=null) {
            timeVaryingReassortmentRates = coalescentDistrInput.get().timeVaryingReassortmentRatesInput.get();
        	isTimeVarying = true;
        }else {
        	reassortmentRate = coalescentDistrInput.get().reassortmentRateInput.get();
        }
        binomialProb = coalescentDistrInput.get().networkIntervalsInput.get().binomialProbInput.get();
        
        redFactor = coalescentDistrInput.get().redFactorInput.get();
        maxHeightRatio = coalescentDistrInput.get().maxHeightRatioInput.get();
        maxHeight = coalescentDistrInput.get().maxHeightInput.get();
    	
        super.initAndValidate();
    }

    @Override
    public double networkProposal() {

    	return resimulate();
    	
    }

    double resimulate() {
    	network.startEditing(this);
        // get the maximum height of the segment tree roots
        double lociHeights = NetworkStatsLogger.getLociMRCA(network);

    	// get all network edges 
        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        // keep only those that coexist at the time of maxHeight
        List<NetworkEdge> startingEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.parentNode.getHeight()>lociHeights)
                .filter(e -> e.childNode.getHeight()<=lociHeights)
               .collect(Collectors.toList());
        

        if (startingEdges.size()==0)
        	return Double.NEGATIVE_INFINITY;
        
       // simulate the rest of the network starting from mxHeight
        double currentTime = lociHeights;
        double timeUntilNextSample = Double.POSITIVE_INFINITY;
        // get the time when the reassortment rates are reduced
        double recChangeTime = Math.min(lociHeights*maxHeightRatio, maxHeight);

        double redFactor = 1.0;
        do {
            // get the current propensities
            int k = startingEdges.size();

            double currentTransformedTime = populationFunction.getIntensity(currentTime);
            double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
            double timeToNextCoal = populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;

            
            double timeToNextReassortment = k>=1 ? 0 : Double.POSITIVE_INFINITY;
            if (isTimeVarying) {
            	if (redFactor==0) {
            		timeToNextReassortment = Double.POSITIVE_INFINITY;
            		
            	}else {
	                double currentTransformedReaTime = timeVaryingReassortmentRates.getIntensity(currentTime);
	                double transformedTimeToNextRea = Randomizer.nextExponential(k * redFactor);
	            	timeToNextReassortment = timeVaryingReassortmentRates.getInverseIntensity(
	            			transformedTimeToNextRea + currentTransformedReaTime) - currentTime;
            	}
            	
            }else {
            	timeToNextReassortment = Randomizer.nextExponential(k*reassortmentRate.getArrayValue() * redFactor);
            }

            // next event time
            double timeUntilNextEvent = Math.min(timeToNextCoal, timeToNextReassortment);
            if ((timeUntilNextEvent+currentTime)>recChangeTime) {
                currentTime = recChangeTime;
                redFactor *= this.redFactor;
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
    
    private void reassort(double reassortmentTime, List<NetworkEdge> extantLineages) {
        NetworkEdge lineage = extantLineages.get(Randomizer.nextInt(extantLineages.size()));

        BitSet hasSegs_left = new BitSet();
        BitSet hasSegs_right = new BitSet();

        for (int segIdx = lineage.hasSegments.nextSetBit(0);
             segIdx != -1; segIdx = lineage.hasSegments.nextSetBit(segIdx+1)) {
        	if (binomialProb==null) {
	            if (Randomizer.nextBoolean()) {
	                hasSegs_left.set(segIdx);
	            } else {
	                hasSegs_right.set(segIdx);
	            }
        	}else {
	            if (Randomizer.nextDouble()>binomialProb.getArrayValue()) {
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
