package coalre.distribution;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.evolution.tree.coalescent.PopulationFunction;

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


    public PopulationFunction populationFunction;
    private Function reassortmentRate;
    public PopulationFunction timeVaryingReassortmentRates;

    public NetworkIntervals intervals;
    
    private boolean isTimeVarying = false;

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

    public double calculateLogP() {
    	logP = 0;
    	// Calculate tree intervals
    	List<NetworkEvent> networkEventList = intervals.getNetworkEventList();

    	NetworkEvent prevEvent = null;

    	for (NetworkEvent event : networkEventList) {
        	if (prevEvent != null)
        		logP += intervalContribution(prevEvent, event);

        	switch (event.type) {
				case COALESCENCE:
					logP += coalesce(event);
					break;

				case SAMPLE:
					break;

				case REASSORTMENT:
					logP += reassortment(event);
					break;
			}

       		if (logP==Double.NEGATIVE_INFINITY)
       			break;

        	prevEvent = event;
        }        
		return logP;
    }
    
	private double reassortment(NetworkEvent event) {
//        lp+=Math.log(reassortmentRate.getArrayValue())
//                + event.segsSortedLeft * Math.log(intervals.getBinomialProb())
//                + (event.segsToSort-event.segsSortedLeft)*Math.log(1-intervals.getBinomialProb())
//                + Math.log(2.0);
		
		double binomval = Math.pow(intervals.getBinomialProb(), event.segsSortedLeft)
				* Math.pow(1-intervals.getBinomialProb(), event.segsToSort-event.segsSortedLeft) 
				+ Math.pow(intervals.getBinomialProb(), event.segsToSort-event.segsSortedLeft)
				* Math.pow(1-intervals.getBinomialProb(), event.segsSortedLeft); 
				        
		if (isTimeVarying)
			return Math.log(timeVaryingReassortmentRates.getPopSize(event.time))
					+ Math.log(binomval);
		else
			return Math.log(reassortmentRate.getArrayValue())
					+ Math.log(binomval);



//        return Math.log(reassortmentRate.getArrayValue())
//                + event.segsSortedLeft * Math.log(intervals.getBinomialProb())
//                + (event.segsToSort-event.segsSortedLeft)*Math.log(1-intervals.getBinomialProb())
//                + Math.log(2.0);
	}

	private double coalesce(NetworkEvent event) {

		return Math.log(1.0/populationFunction.getPopSize(event.time));
	}

	private double intervalContribution(NetworkEvent prevEvent, NetworkEvent nextEvent) {

        double result = 0.0;

//        System.out.println(timeVaryingReassortmentRates.getIntegral(prevEvent.time, nextEvent.time));
//        System.out.println(nextEvent.time - prevEvent.time);
		if (isTimeVarying)
	        result += -prevEvent.totalReassortmentObsProb 
	        	* timeVaryingReassortmentRates.getIntegral(prevEvent.time, nextEvent.time);
		else
	        result += -reassortmentRate.getArrayValue() * prevEvent.totalReassortmentObsProb
	                * (nextEvent.time - prevEvent.time);

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
