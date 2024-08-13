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
public class PiecewiseConstantReassortmentRateScalers extends PopulationFunction.Abstract {
    final public Input<RealParameter> InfectedToRhoInput = new Input<>("InfectedToRho",
            "the value that maps the number of infected or the Ne to the reassortment rate ");
    final public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts",
            "When to switch between elements of Ne", Input.Validate.REQUIRED);

    
    final public Input<Spline> splineInput = new Input<>("spline",
            "Spline to use for the population function", Input.Validate.REQUIRED);

    boolean NesKnown = false;

    Spline spline;
    RealParameter InfectedToRho;
    RealParameter rateShifts;
    

    double rateRatio;
    
    boolean ratioIsSpline = false;

    @Override
    public void initAndValidate() {
        rateShifts = rateShiftsInput.get();
        InfectedToRho = InfectedToRhoInput.get();
        // ensure that the reassortment rates is of dimension rateShifts
		if (InfectedToRho.getDimension() != rateShifts.getDimension()) {
			InfectedToRho.setDimension(rateShifts.getDimension());
		}

        spline = splineInput.get();
    }


    @Override
    public List<String> getParameterIds() {
        return null;
    }

    @Override
    public double getPopSize(double t) {
        // check which time t is, if it is larger than the last time, return the last Ne
        int interval = spline.gridPoints-1;
        for (int i = 0; i < spline.gridPoints; i++){
            if (t < spline.time[i]){
                interval = i-1;
                break;
            }
        }
        
        for (int i = 1; i < rateShifts.getDimension(); i++) {
			if (t >= rateShifts.getValue(i-1) && t < rateShifts.getValue(i)) {
				return spline.I[interval] * InfectedToRho.getValue(i-1);
            }
        }
        
        return spline.I[interval] * InfectedToRho.getValue(rateShifts.getDimension() - 1);
        
			

    }

    public double getIntegral(double from, double to) {
        // compute the integral of Ne's between from an to
        double NeIntegral = 0;
        int intervalFrom = spline.gridPoints-1;
        for (int i = 0; i < spline.gridPoints; i++) {
            // get the first time larger than from
            if (from < spline.time[i]) {
                intervalFrom = i-1;
                break;
            }
        }
        for (int i = intervalFrom; i < spline.gridPoints; i++) {
        	double rateRatio=0.0;    
        	boolean found = false;
            for (int j = 1; j < rateShifts.getDimension(); j++) {
    			if (spline.time[i] >= rateShifts.getValue(j-1) && spline.time[i] < rateShifts.getValue(j)) {
    				rateRatio = InfectedToRho.getValue(j-1);
    				found = true;
                }
            }
			if (!found)
				rateRatio = InfectedToRho.getValue(rateShifts.getDimension() - 1);            
        	
        	
            double rate = spline.I[i] * rateRatio;
			
            // if i==intervalFrom, we have to start compute the diff from there
            if (spline.time[i+1] > to) {
                if (i == intervalFrom) {
                    NeIntegral += rate * (to - from);
                }else{
                    NeIntegral += rate * (to - spline.time[i]);
                }
                break;
            }else if (i == intervalFrom) {
                NeIntegral += rate * (spline.time[i+1] - from);
            }else{
                NeIntegral += rate * (spline.time[i + 1] - spline.time[i]);
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
		throw new IllegalArgumentException("Not implemented");
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