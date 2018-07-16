package coalre.distribution;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;

import java.util.BitSet;
import java.util.List;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a reassortment network using under" +
        " the framework of Mueller (2018).")
public class CoalescentWithReassortment extends NetworkDistribution {
	
	public Input<RealParameter> reassortmentRateInput = new Input<>(
	        "reassortmentRate",
            "reassortment rate (per lineage per unit time)",
            Input.Validate.REQUIRED);

	public Input<PopulationFunction> populationFunctionInput = new Input<>(
	        "populationModel",
            "Population model.",
            Input.Validate.REQUIRED);

	public Input<Boolean> simpleReassortmentOnlyInput = new Input<>(
	        "simpleReassortmentOnly",
            "Disable generation of reassortment events on lineages not ancestral to all segments.",
            false);
	
    private PopulationFunction populationFunction;
    private RealParameter reassortmentRate;
    private boolean simpleReassortmentOnly;

    @Override
    public void initAndValidate(){
        populationFunction = populationFunctionInput.get();
        reassortmentRate = reassortmentRateInput.get();
        simpleReassortmentOnly = simpleReassortmentOnlyInput.get();
    }

    public double calculateLogP() {
    	logP = 0;

    	// Calculate tree intervals
    	List<NetworkEvent> networkEventList = networkIntervalsInput.get().getNetworkEventList();

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

    	BitSet segments = event.node.getChildEdges().get(0).hasSegments;

    	if (simpleReassortmentOnly && segments.cardinality()<networkIntervalsInput.get().getSegmentCount())
    	    return Double.NEGATIVE_INFINITY;

        return Math.log(reassortmentRate.getValue()) +
                Math.log(1.0/segments.cardinality()) +
                Math.log(1.0/(segments.cardinality()-1)) +
                (segments.cardinality()-2)*Math.log(0.5);
	}

	private double coalesce(NetworkEvent event) {

		return Math.log(1.0/populationFunction.getPopSize(event.time));
	}

	private double intervalContribution(NetworkEvent prevEvent, NetworkEvent nextEvent) {

        double result = 0.0;

        result += -reassortmentRate.getValue()*prevEvent.totalReassortmentObsProb
                * (nextEvent.time-prevEvent.time);

		
		result += -0.5*prevEvent.lineages*(prevEvent.lineages-1)
                * populationFunction.getIntegral(prevEvent.time, nextEvent.time);
		
		return result;
	}   

    
}
