package coalre.distribution;

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
            "reassortment rate (per lineage per unit time)",
            Input.Validate.REQUIRED);

	public Input<PopulationFunction> populationFunctionInput = new Input<>(
	        "populationModel",
            "Population model.",
            Input.Validate.REQUIRED);

    private PopulationFunction populationFunction;
    private Function reassortmentRate;
    private NetworkIntervals intervals;

    @Override
    public void initAndValidate(){
        populationFunction = populationFunctionInput.get();
        reassortmentRate = reassortmentRateInput.get();
        intervals = networkIntervalsInput.get();
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

        return Math.log(reassortmentRate.getArrayValue())
                + event.segsSortedLeft * Math.log(intervals.getBinomialProb())
                + (event.segsToSort-event.segsSortedLeft)*Math.log(1-intervals.getBinomialProb())
                + Math.log(2.0);
	}

	private double coalesce(NetworkEvent event) {

		return Math.log(1.0/populationFunction.getPopSize(event.time));
	}

	private double intervalContribution(NetworkEvent prevEvent, NetworkEvent nextEvent) {

        double result = 0.0;

        result += -reassortmentRate.getArrayValue() * prevEvent.totalReassortmentObsProb
                * (nextEvent.time - prevEvent.time);

		result += -0.5*prevEvent.lineages*(prevEvent.lineages-1)
                * populationFunction.getIntegral(prevEvent.time, nextEvent.time);
		
		return result;
	}
}
