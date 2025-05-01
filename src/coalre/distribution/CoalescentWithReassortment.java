package coalre.distribution;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import coalre.statistics.NetworkStatsLogger;

import java.util.List;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a reassortment network using under" +
        " the framework of Mueller (2018).")
public class CoalescentWithReassortment extends NetworkDistribution {
	
	public Input<Function> reassortmentRateInput = new Input<>(
	        "reassortmentRate",
            "reassortment rate (per lineage per unit time)");

	public Input<PopulationFunction> populationFunctionInput = new Input<>(
	        "populationModel",
            "Population model.",
            Input.Validate.REQUIRED);
	
	public Input<PopulationFunction> timeVaryingReassortmentRatesInput = new Input<>(
	        "timeVaryingReassortmentRates",
            "reassortment rates that vary over time",
            Input.Validate.XOR, reassortmentRateInput);

	public Input<Double> maxHeightRatioInput = new Input<>(
			"maxHeightRatio",
			"if specified, above the ratio, only coalescent events are allowed.", Double.POSITIVE_INFINITY);

	public Input<Double> redFactorInput = new Input<>(
			"redFactor",
			"by how much the recombination rate should be reduced after reaching the maxHeightRatio.", 0.1);



	public PopulationFunction populationFunction;
    private Function reassortmentRate;
    public PopulationFunction timeVaryingReassortmentRates;

    public NetworkIntervals intervals;
    
    private boolean isTimeVarying = false;

	public double redFactor;


	@Override
    public void initAndValidate(){
        populationFunction = populationFunctionInput.get();
        intervals = networkIntervalsInput.get();
        if (reassortmentRateInput.get()!=null) {
        	reassortmentRate = reassortmentRateInput.get();
        }else {
        	isTimeVarying = true;
        	timeVaryingReassortmentRates = timeVaryingReassortmentRatesInput.get();
        }
    }
	int ii = 0;
    public double calculateLogP() {
    	logP = 0;
    	// Calculate tree intervals
    	List<NetworkEvent> networkEventList = intervals.getNetworkEventList();

		// get the mrca of all loci trees
		double lociMRCA = maxHeightRatioInput.get()<Double.POSITIVE_INFINITY ?
				NetworkStatsLogger.getLociMRCA(networkIntervalsInput.get().networkInput.get()) : Double.POSITIVE_INFINITY;
		
//		System.out.println(networkIntervalsInput.get().networkInput.get());
//		System.out.println("lociMRCA: " + NetworkStatsLogger.getLociMRCA(networkIntervalsInput.get().networkInput.get()));
		
		NetworkEvent prevEvent = null;

    	for (NetworkEvent event : networkEventList) {
        	if (prevEvent != null)
        		logP += intervalContribution(prevEvent, event, lociMRCA);

        	switch (event.type) {
				case COALESCENCE:
					logP += coalesce(event);
					break;

				case SAMPLE:
					break;

				case REASSORTMENT:
					logP += reassortment(event, lociMRCA);
					break;
			}

        	

       		if (logP==Double.NEGATIVE_INFINITY)
       			break;

        	prevEvent = event;
        }
//    	System.out.println("logP: " + logP);
//    	ii++;
//    	if (ii>100) {
//    		System.exit(0);
//    	}
//    	if (logP>0) {
//    		prevEvent = null;
//
//        	for (NetworkEvent event : networkEventList) {
//            	if (prevEvent != null)
//            		System.out.println(intervalContribution(prevEvent, event, lociMRCA));
//            	switch (event.type) {
//    				case COALESCENCE:
//    					System.out.println("c" +coalesce(event));
//    					break;
//
//    				case SAMPLE:
//    					break;
//
//    				case REASSORTMENT:
//    					System.out.println("r" +reassortment(event, lociMRCA));
//    					break;
//            	}
//            	prevEvent = event;
//        	}
//
//    		
//    		
//    		System.out.println(networkIntervalsInput.get().networkInput.get());
//    		System.exit(0);
//    	}
//		System.out.println(networkIntervalsInput.get().networkInput.get());
//		System.exit(0);
//		System.out.println(logP);
//		if (logP == Double.POSITIVE_INFINITY)
//			System.exit(0);

		return logP;
    }
    
	private double reassortment(NetworkEvent event, double lociMRCA) {
		
		double binomval = Math.pow(intervals.getBinomialProb(), event.segsSortedLeft)
				* Math.pow(1-intervals.getBinomialProb(), event.segsToSort-event.segsSortedLeft) 
				+ Math.pow(intervals.getBinomialProb(), event.segsToSort-event.segsSortedLeft)
				* Math.pow(1-intervals.getBinomialProb(), event.segsSortedLeft);
		
		if (event.time<=(lociMRCA*maxHeightRatioInput.get())) {
			if (isTimeVarying) {
				return Math.log(timeVaryingReassortmentRates.getPopSize(event.time))
						+ Math.log(binomval);
				
			}else
				return Math.log(reassortmentRate.getArrayValue())
						+ Math.log(binomval);
		}else{
			if (isTimeVarying)
				return Math.log(redFactor*timeVaryingReassortmentRates.getPopSize(event.time))
						+ Math.log(binomval);
			else
				return Math.log(redFactor*reassortmentRate.getArrayValue())
						+ Math.log(binomval);
		}	
	}

	private double coalesce(NetworkEvent event) {
		return Math.log(1.0/populationFunction.getPopSize(event.time));
	}

	private double intervalContribution(NetworkEvent prevEvent, NetworkEvent nextEvent, double lociMRCA) {

        double result = 0.0;
        

		if (nextEvent.time<(lociMRCA*maxHeightRatioInput.get())) {
			if (isTimeVarying) {
				result += -prevEvent.totalReassortmentObsProb
						* timeVaryingReassortmentRates.getIntegral(prevEvent.time, nextEvent.time);
			}else {
				result += -reassortmentRate.getArrayValue() * prevEvent.totalReassortmentObsProb
						* (nextEvent.time - prevEvent.time);

			}
		}else if (prevEvent.time<(lociMRCA*maxHeightRatioInput.get())) {
			if (isTimeVarying) {
				result += -prevEvent.totalReassortmentObsProb
						* timeVaryingReassortmentRates.getIntegral(prevEvent.time, lociMRCA*maxHeightRatioInput.get());
			}else {
				result += -reassortmentRate.getArrayValue() * prevEvent.totalReassortmentObsProb
						* (lociMRCA*maxHeightRatioInput.get() - prevEvent.time);
			}

			if (isTimeVarying) {
				result += -redFactor*prevEvent.totalReassortmentObsProb
						* timeVaryingReassortmentRates.getIntegral(lociMRCA*maxHeightRatioInput.get(), nextEvent.time);
			}else {
				result += -redFactor*reassortmentRate.getArrayValue() * prevEvent.totalReassortmentObsProb
						* (nextEvent.time - prevEvent.time-lociMRCA*maxHeightRatioInput.get());
			}



		}else{
			if (isTimeVarying)
				result += -prevEvent.totalReassortmentObsProb
						* redFactor* timeVaryingReassortmentRates.getIntegral(prevEvent.time, nextEvent.time);
			else
				result += -redFactor * reassortmentRate.getArrayValue() * prevEvent.totalReassortmentObsProb
						* (nextEvent.time - prevEvent.time);
		}
		result += -0.5*prevEvent.lineages*(prevEvent.lineages-1)
                * populationFunction.getIntegral(prevEvent.time, nextEvent.time);
		
		return result;
	}
	
    @Override
    protected boolean requiresRecalculation() {    	
    	if (isTimeVarying) {
	    	if (((CalculationNode) timeVaryingReassortmentRates).isDirtyCalculation())
	    		return true;
    	}else{
	    	if (((CalculationNode) reassortmentRate).isDirtyCalculation())
	    		return true;
    	}
    	if (((CalculationNode) populationFunction).isDirtyCalculation())
    		return true;
    	
        return super.requiresRecalculation();
    }


    

}
