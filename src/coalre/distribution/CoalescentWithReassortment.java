package coalre.distribution;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import coalre.network.NetworkNode;
import coalre.statistics.NetworkStatsLogger;

import java.util.Arrays;
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
	
	public Input<Double> maxHeightInput = new Input<>(
			"maxHeight",
			"if specified, above the height, only coalescent events are allowed.", Double.POSITIVE_INFINITY);

	public Input<Double> redFactorInput = new Input<>(
			"redFactor",
			"by how much the recombination rate should be reduced after reaching the maxHeightRatio.", 0.1);



	public PopulationFunction populationFunction;
    private Function reassortmentRate;
    public PopulationFunction timeVaryingReassortmentRates;

    public NetworkIntervals intervals;
    
    private boolean isTimeVarying = false;

	public double redFactor;
	double[][] logBinomval;
		
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
        
        logBinomval = new double[intervals.networkInput.get().getSegmentCount()+1][intervals.networkInput.get().getSegmentCount()+1];
		for (int i = 0; i < intervals.networkInput.get().getSegmentCount()+1; i++) {
			
			for (int j = 0; j <= i; j++) {
				if (i==0 && j==0) {
					logBinomval[i][j] = Math.log(0.1*intervals.emptyObsProb);
				}else if (i==j || j==0) {
					logBinomval[i][j] = Math.log((i+1)*intervals.emptyObsProb*(Math.pow(intervals.getBinomialProb(), j)
							* Math.pow(1-intervals.getBinomialProb(), i-j) 
							+ Math.pow(intervals.getBinomialProb(), i-j)
							* Math.pow(1-intervals.getBinomialProb(), j)));
//
//				}else if (j==0) {
//					logBinomval[i][j] = Math.pow(intervals.getBinomialProb(), i);
				}else {
					logBinomval[i][j] = Math.log(Math.pow(intervals.getBinomialProb(), j)
							* Math.pow(1-intervals.getBinomialProb(), i-j) 
							+ Math.pow(intervals.getBinomialProb(), i-j)
							* Math.pow(1-intervals.getBinomialProb(), j));
				}
			}
//			System.out.println(Arrays.toString(logBinomval[i]));
		}
//        System.exit(0);
    }

	public double calculateLogP() {		
    	logP = 0;
    	// Calculate tree intervals
    	intervals.update();
		// get the mrca of all loci trees
		double oriLociMRCA = NetworkStatsLogger.getLociMRCA(networkIntervalsInput.get().networkInput.get());
		
		double lociMRCA = oriLociMRCA;
		// get the second highest segment root
		if (maxHeightRatioInput.get()<Double.POSITIVE_INFINITY || maxHeightInput.get()<Double.POSITIVE_INFINITY)
			lociMRCA = NetworkStatsLogger.getSecondHighestLociMRCA(networkIntervalsInput.get().networkInput.get());
		
		double maxHeight = Math.min(maxHeightInput.get(),lociMRCA * maxHeightRatioInput.get());

		NetworkEvent prevEvent = null;
		
		double logPint = 0;
		double logPcoal = 0.;
		double logPreassortment = 0.;

    	for (NetworkEvent event : intervals.networkEventList) {
        	if (prevEvent != null) {
        		logPint += intervalContribution(prevEvent, event, maxHeight, oriLociMRCA);
//        		logP += logPint;
        	}
        	
        	switch (event.type) {
				case COALESCENCE:
					logPcoal += coalesce(event);
//					logP += logPcoal;
					break;

				case SAMPLE:
					break;

				case REASSORTMENT:
					logPreassortment += reassortment(event, maxHeight);
//					logP += logPreassortment;
					break;
			}


       		if (logP==Double.NEGATIVE_INFINITY)
       			break;

        	prevEvent = event;
        }
    	logP = logPint + logPcoal + logPreassortment;
    	if (logP>0) {
    		System.err.println(logP + " " + logPint + " " + logPcoal + " " + logPreassortment);
    		System.exit(0);
    	}

		return logP;
    }
    
	private double reassortment(NetworkEvent event, double maxHeight) {
		
		double logBinomval;
		if (intervals.getBinomialProb()!=0.5) {
			logBinomval = Math.log(Math.pow(intervals.getBinomialProb(), event.segsSortedLeft)
				* Math.pow(1-intervals.getBinomialProb(), event.segsToSort-event.segsSortedLeft) 
				+ Math.pow(intervals.getBinomialProb(), event.segsToSort-event.segsSortedLeft)
				* Math.pow(1-intervals.getBinomialProb(), event.segsSortedLeft));
		}else {
			try {
				logBinomval = this.logBinomval[event.segsToSort][event.segsSortedLeft];
			} catch (ArrayIndexOutOfBoundsException e) {
				System.err.println(intervals.networkInput.get());
				System.err.println(event.segsToSort + " " + event.segsSortedLeft);
				for (NetworkNode n: intervals.networkInput.get().getNodes()) {
					if (n.getHeight() == event.time) {
						System.err.println(n.getChildCount());
						System.err.println(n.getParentCount());
					}
				}
				System.err.println("empty edge remaining, should not happen, event at time "+ event.time);
				return Double.NEGATIVE_INFINITY;
//				logBinomval=1;
			}
		}
		
		if (event.time<=maxHeight) {
			if (isTimeVarying) {
				return Math.log(timeVaryingReassortmentRates.getPopSize(event.time))
						+ logBinomval;
				
			}else
				return Math.log(reassortmentRate.getArrayValue())
						+ logBinomval;
		}else{
			if (isTimeVarying)
				return Math.log(redFactor*timeVaryingReassortmentRates.getPopSize(event.time))
						+ logBinomval;
			else
				return Math.log(redFactor*reassortmentRate.getArrayValue())
						+ logBinomval;
		}	
	}

	private double coalesce(NetworkEvent event) {
		return Math.log(1.0/populationFunction.getPopSize(event.time));
	}

	private double intervalContribution(NetworkEvent prevEvent, NetworkEvent nextEvent, double maxHeight, double oriLociMRCA) {

        double result = 0.0;
               
		if (nextEvent.time<maxHeight) {
			if (isTimeVarying) {
				result += -prevEvent.totalReassortmentObsProb
						* timeVaryingReassortmentRates.getIntegral(prevEvent.time, nextEvent.time);
			}else {
				result += -reassortmentRate.getArrayValue() * prevEvent.totalReassortmentObsProb
						* (nextEvent.time - prevEvent.time);

			}
		}else if (prevEvent.time<maxHeight) {
			if (isTimeVarying) {
				result += -prevEvent.totalReassortmentObsProb
						* timeVaryingReassortmentRates.getIntegral(prevEvent.time, maxHeight);
			}else {
				result += -reassortmentRate.getArrayValue() * prevEvent.totalReassortmentObsProb
						* (maxHeight - prevEvent.time);
			}

			if (isTimeVarying) {
				result += -redFactor*prevEvent.totalReassortmentObsProb
						* timeVaryingReassortmentRates.getIntegral(maxHeight, nextEvent.time);
			}else {
				result += -redFactor*reassortmentRate.getArrayValue() * prevEvent.totalReassortmentObsProb
						* (nextEvent.time - prevEvent.time-maxHeight);
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
		
		if (nextEvent.time > oriLociMRCA && prevEvent.lineages == 1) {
//			System.out.println(intervals.networkInput.get());
//			System.err.println("Warning");
			// ensures that root is actually the root, this was previous done by the network terminates at tMRCA class
			// however, this is a more efficient way of keeping track of that
			return Double.NEGATIVE_INFINITY;
		}
		
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
