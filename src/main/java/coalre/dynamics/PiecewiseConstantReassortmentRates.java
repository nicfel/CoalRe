package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;

import java.io.PrintStream;
import java.util.List;


/**
 * @author Nicola F. Mueller
 */
@Description("Computes time varying recombination rates as a population"+
        " function from spline interpolation of the number of infected over time")
public class PiecewiseConstantReassortmentRates extends PopulationFunction.Abstract {
    final public Input<RealParameter> reassortmentRateInput = new Input<>("reassortmentRate",
            "the value of the reassortmentRate");
    final public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts",
            "When to switch between elements of Ne", Input.Validate.REQUIRED);


    RealParameter reassortmentRate;
    RealParameter rateShifts;


    double rateRatio;
    
    boolean ratioIsSpline = false;

    @Override
    public void initAndValidate() {
        rateShifts = rateShiftsInput.get();
        reassortmentRate = reassortmentRateInput.get();
        // ensure that the reassortment rates is of dimension rateShifts
		if (reassortmentRate.getDimension() != rateShifts.getDimension()) {
			reassortmentRate.setDimension(rateShifts.getDimension());
		}
    }


    @Override
    public List<String> getParameterIds() {
        return null;
    }

    @Override
    public double getPopSize(double t) {
    	// check which time t is, if it is larger than the last time, return the last Ne
        for (int i = 1; i < rateShifts.getDimension(); i++) {
			if (t >= rateShifts.getValue(i-1) && t < rateShifts.getValue(i)) {
				return Math.exp(reassortmentRate.getValue(i-1));
            }
        }
        return Math.exp(reassortmentRate.getValue(rateShifts.getDimension() - 1));
        
    }

    public double getIntegral(double from, double to) {

        // compute the integral of Ne's between from an to
        double NeIntegral = 0;
        double fromTime = 0, toTime=0;
		for (int i = 1; i < rateShifts.getDimension(); i++) {
			// check if from is in this interval
			if (from >= rateShifts.getValue(i - 1) && from < rateShifts.getValue(i)) {
				fromTime = from;
			}else if (from < rateShifts.getValue(i - 1)) {
				fromTime = rateShifts.getValue(i - 1);
            }else{
            	fromTime = rateShifts.getValue(i);
            }
			
			// check if to in in this interval
			if (to >= rateShifts.getValue(i - 1) && to < rateShifts.getValue(i)) {
				toTime = to;
			} else if (to < rateShifts.getValue(i - 1)) {
				break;
			} else {
				toTime = rateShifts.getValue(i);
			}
			NeIntegral += (toTime - fromTime) * Math.exp(reassortmentRate.getValue(i - 1));
		}       

		// add the scenario where to or to and from are larger than the last time
		if (to >= rateShifts.getValue(rateShifts.getDimension() - 1)) {
			if (from >= rateShifts.getValue(rateShifts.getDimension() - 1)) {
				NeIntegral += (to - from) * Math.exp(reassortmentRate.getValue(rateShifts.getDimension() - 1));
			} else {
				NeIntegral += (to - rateShifts.getValue(rateShifts.getDimension() - 1))
						* Math.exp(reassortmentRate.getValue(rateShifts.getDimension() - 1));
			}
		}
		
        return NeIntegral;
    }

    @Override
    public double getIntensity(double v){
        return getIntegral(0, v);
    }

    @Override
    public double getInverseIntensity(double v) {
    	
        // compute the integral of Ne's between from an to
		for (int i = 1; i < rateShifts.getDimension(); i++) {
			v -= (rateShifts.getValue(i) - rateShifts.getValue(i - 1)) * Math.exp(reassortmentRate.getValue(i - 1));
			if (v < 0) {
				v += (rateShifts.getValue(i) - rateShifts.getValue(i - 1)) * Math.exp(reassortmentRate.getValue(i - 1));
				// solve for the final time
				return rateShifts.getValue(i - 1) + v / Math.exp(reassortmentRate.getValue(i - 1));
			}			
		}    		

		return rateShifts.getValue(rateShifts.getDimension() - 1)
				+ v / Math.exp(reassortmentRate.getValue(rateShifts.getDimension() - 1));
    }

    @Override
    public boolean requiresRecalculation() {
        return true;
    }

    @Override
    public void store() {
        super.store();
    }

    @Override
    public void restore() {
        super.restore();
    }



}