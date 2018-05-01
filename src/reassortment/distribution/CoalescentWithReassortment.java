package reassortment.distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


import beast.core.CalculationNode;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.coalescent.IntervalType;
import beast.util.Randomizer;
import reassortment.network.NetworkIntervalType;
import reassortment.network.Networknode;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a beast.tree using under the framework of Mueller (2017).")
public class CoalescentWithReassortment extends NetworkDistribution {
	
	public Input<RealParameter> rRateInput = new Input<>("reassortmentRate", "reassortment rate", Input.Validate.REQUIRED);
	public Input<RealParameter> coalescentRateInput = new Input<>("coalescentRate", "coalescent rate", Input.Validate.REQUIRED);
	

	int failedParticles;
    
    // Set up for lineage state probabilities
    private ArrayList<Integer> activeLineages;
    private ArrayList<Boolean[]> activeSegments;
    
    private double coalescentRate;
    private double reassortmentRate;
           
    @Override
    public void initAndValidate(){    	
    }
    
//    int count = 0;
    public double calculateLogP() {
//    	System.out.println("calc");
    	logP = 0;
    	// newly calculate tree intervals
    	networkIntervalsInput.get().calculateIntervals();
    	double nextEventTime = 0.0;
    	int networkInterval = 0;
    	
    	boolean first = true;
    	
    	activeLineages = new ArrayList<>();
    	activeSegments = new ArrayList<>();
    	
    	coalescentRate = coalescentRateInput.get().getValue(0);
    	reassortmentRate = rRateInput.get().getValue(0);
    	
    	
    	
        do {       
        	nextEventTime = networkIntervalsInput.get().getInterval(networkInterval);
        	
        	if (nextEventTime==Double.POSITIVE_INFINITY)
        		break;
       	
        	if (nextEventTime > 0) {
        		logP += intervalContribution(nextEventTime);
        	}
       	
        	if (networkIntervalsInput.get().getNetworkIntervalType(networkInterval) == NetworkIntervalType.COALESCENT) {
//        		System.out.println("coal");
        		logP += coalesce(networkInterval);
        	}
       		
       		if (networkIntervalsInput.get().getNetworkIntervalType(networkInterval) == NetworkIntervalType.SAMPLE) { 	
//        		System.out.println("sam");
      			sample(networkInterval, first);
       			first = false;
       		}	 
       		
       		if (networkIntervalsInput.get().getNetworkIntervalType(networkInterval) == NetworkIntervalType.REASSORTMENTSTART) { 
//        		System.out.println("REASSORTMENTSTART");
        		logP += reassortmentstart(networkInterval);
       			first = false;
       		}	 
       		
       		if (networkIntervalsInput.get().getNetworkIntervalType(networkInterval) == NetworkIntervalType.REASSORTMENTEND) { 	
//        		System.out.println("REASSORTMENTEND");
      			reassortmentend(networkInterval);
       			first = false;
       		}	    		
       		
       		
       		networkInterval++;
       		if (logP==Double.NEGATIVE_INFINITY)
       			break;
        }while(nextEventTime <= Double.POSITIVE_INFINITY);
        
//        count++;
        
//        System.out.println(logP);
//        if (count == 2)
//        	System.exit(0);
		return logP;  	
    }
    
	private double reassortmentstart(int networkInterval) {
		List<Networknode> reassLines = networkIntervalsInput.get().getLineagesRemoved(networkInterval);

    	if (reassLines.size() != 1) {
			System.err.println("Unsupported number of incoming lineages at a reassortment event");
			System.exit(0);
		}
    	
    	final int daughterIndex = activeLineages.indexOf(reassLines.get(0).getNr());
    	
		if (daughterIndex == -1 ) {
			System.out.println(reassLines.get(0).getNr() + " " + activeLineages);
			System.out.println("daughter lineage at reassortment event not found");
			return Double.NaN;
		}


    	Boolean[] segsBefore = reassLines.get(0).getHasSegments();
    	
    	int c_before = 0;
    	for (int i = 0; i < segsBefore.length; i++){
    		if (segsBefore[i]==true) c_before++;
    	}
    	
    	
        activeLineages.add(reassLines.get(0).getParent().getNr());  
        
        activeSegments.add(reassLines.get(0).getParent().getHasSegments());
        
        activeLineages.remove(daughterIndex);
        activeSegments.remove(daughterIndex);


		return Math.log(reassortmentRate*2*Math.pow(0.5, c_before));
		
	}


	private void reassortmentend(int networkInterval) {
		List<Networknode> reassLines = networkIntervalsInput.get().getLineagesAdded(networkInterval);
		
    	if (reassLines.size() != 1) {
			System.err.println("Unsupported number of incoming lineages at a second part of a reassortment event");
			System.exit(0);
    	}
		
		activeLineages.add(reassLines.get(0).getNr());
		activeSegments.add(reassLines.get(0).getHasSegments());
	}


	private void sample(int networkInterval, boolean first) {
		List<Networknode> incomingLines = networkIntervalsInput.get().getLineagesAdded(networkInterval);	
		
		for (Networknode l : incomingLines) {
			activeLineages.add(l.getNr());
			activeSegments.add(l.getHasSegments());
		}
	}

	private double coalesce(int networkInterval) {
		List<Networknode> coalLines = networkIntervalsInput.get().getLineagesRemoved(networkInterval);
		
    	if (coalLines.size() > 2) {
			System.err.println("Unsupported coalescent at non-binary node");
			System.exit(0);
		}
    	if (coalLines.size() < 2) {
    		System.out.println();
    		System.out.println("WARNING: Less than two lineages found at coalescent event!");
    		System.out.println();
    		return Double.NaN;
		}
		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println(coalLines.get(0).getNr() + " " + coalLines.get(1).getNr() + " " + activeLineages);
			System.out.println("daughter lineages at coalescent event not found");
			return Double.NaN;
		}
		
        activeLineages.add(coalLines.get(0).getParent().getNr());    
        
        activeSegments.add(coalLines.get(0).getParent().getHasSegments());
        
		activeLineages.remove(Math.max(daughterIndex1, daughterIndex2));
		activeLineages.remove(Math.min(daughterIndex1, daughterIndex2));
		
		activeSegments.remove(Math.max(daughterIndex1, daughterIndex2));
		activeSegments.remove(Math.min(daughterIndex1, daughterIndex2));

		
		return Math.log(coalescentRate);
	}

	private double intervalContribution(double nextEventTime) {
		// for each active lineage, calculate the reassortment interval contribution
		double reassCont = 0.0;
		// TODO make this thing more efficient
		for (int i = 0; i < activeSegments.size(); i++){
	    	int nr_lins = 0;
	    	for (int j = 0; j < activeSegments.get(i).length; j++)
	    		if (activeSegments.get(i)[j]==true) nr_lins++;
	    	
	    	// Precompute powers
	    	reassCont -= (1-2*Math.pow(0.5, nr_lins));

		}
		reassCont *= reassortmentRate;
		
		
		int overlap = 0;
		
		// TODO keep directly track of which network nodes overlap
		for (int i = 0; i < activeSegments.size(); i++){
			for (int j = i+1; j < activeSegments.size(); j++){
				for (int k = 0; k < activeSegments.get(i).length; k++){
					if (activeSegments.get(i)[k] && activeSegments.get(i)[k]){
						overlap++;
						break; // break out of the loop if there is at least one segment overlapping						
					}// if
				}// k					
			}// j		
		}// i
		
		double coalCont = -overlap*coalescentRate;
		
		return nextEventTime*(reassCont+coalCont);
	}   

    
}
