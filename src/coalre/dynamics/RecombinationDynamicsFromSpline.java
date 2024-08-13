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
public class RecombinationDynamicsFromSpline extends PopulationFunction.Abstract implements Loggable {
    final public Input<RealParameter> InfectedToRhoInput = new Input<>("InfectedToRho",
            "the value that maps the number of infected or the Ne to the reassortment rate ");
    final public Input<Spline> InfectedToRhoSplineInput = new Input<>("infectedToRhoSpline",
            "the value that maps the number of infected or the Ne to the reassortment rate ", Input.Validate.XOR, InfectedToRhoInput);
       
    
    final public Input<Spline> splineInput = new Input<>("spline",
            "Spline to use for the population function");

    boolean NesKnown = false;

    Spline spline;
    Spline infectedToRhoSpline;    
    RealParameter InfectedToRho;
    boolean hasNeSpline = false;
    

    double rateRatio;    
    boolean ratioIsSpline = false;

    @Override
    public void initAndValidate() {
    	if (InfectedToRhoSplineInput.get() != null) {
    		infectedToRhoSpline = InfectedToRhoSplineInput.get();
    		ratioIsSpline = true;
    	} else {
    		InfectedToRho = InfectedToRhoInput.get();
    	}
    	if (splineInput.get() != null) {
			spline = splineInput.get();
			hasNeSpline = true;
    	}
    }


    @Override
    public List<String> getParameterIds() {
        return null;
    }

    @Override
    public double getPopSize(double t) {
		if (hasNeSpline) {

	        if (!spline.update())
	            return Double.NaN;
			if (ratioIsSpline && !infectedToRhoSpline.update())
				return Double.NaN;
	
	        // check which time t is, if it is larger than the last time, return the last Ne
	        int interval = spline.gridPoints-1;
	        for (int i = 0; i < spline.gridPoints; i++){
	            if (t < spline.time[i]){
	                interval = i-1;
	                break;
	            }
	        }
			if (ratioIsSpline)
				return spline.I[interval] * infectedToRhoSpline.I[interval];			
			else
				return spline.I[interval] * InfectedToRho.getValue();
	    }else {
	    	if (!infectedToRhoSpline.update())
				return Double.NaN;

	        // check which time t is, if it is larger than the last time, return the last Ne
	        int interval = infectedToRhoSpline.gridPoints-1;
	        for (int i = 0; i < infectedToRhoSpline.gridPoints; i++){
	            if (t < infectedToRhoSpline.time[i]){
	                interval = i-1;
	                break;
	            }
	        }

			return infectedToRhoSpline.I[interval];	    	
	    }

    }

    public double getIntegral(double from, double to) {
//    	System.out.println(from + " " + to);
    	if (from<0 && from > -1e-14) {
    		from=0;
    	}
    	if (hasNeSpline) {
	        if (!spline.update())
	            return Double.NaN;
			if (ratioIsSpline && !infectedToRhoSpline.update())
				return Double.NaN;
	
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
				double rate = spline.I[i];
				
				if (ratioIsSpline)
					rate *= infectedToRhoSpline.I[i];
				else
					rate *= InfectedToRho.getValue();
				
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
    	}else {

			if (!infectedToRhoSpline.update())
				return Double.NaN;
	
	        // compute the integral of Ne's between from an to
	        double NeIntegral = 0;
	        int intervalFrom = infectedToRhoSpline.gridPoints-1;
	        for (int i = 0; i < infectedToRhoSpline.gridPoints; i++) {
	            // get the first time larger than from
	            if (from < infectedToRhoSpline.time[i]) {
	                intervalFrom = i-1;
	                break;
	            }
	        }
	        for (int i = intervalFrom; i < infectedToRhoSpline.gridPoints; i++) {
				double rate = infectedToRhoSpline.I[i];
				
	            // if i==intervalFrom, we have to start compute the diff from there
	            if (infectedToRhoSpline.time[i+1] > to) {
	                if (i == intervalFrom) {
	                    NeIntegral += rate * (to - from);
	                }else{
	                    NeIntegral += rate * (to - infectedToRhoSpline.time[i]);
	                }
	                break;
	            }else if (i == intervalFrom) {
	                NeIntegral += rate * (infectedToRhoSpline.time[i+1] - from);
	            }else{
	                NeIntegral += rate * (infectedToRhoSpline.time[i + 1] - infectedToRhoSpline.time[i]);
	            }
	        }
	        return NeIntegral;
    		
    	}
    }

    @Override
    public double getIntensity(double v){
        return getIntegral(0, v);
    }

    @Override
    public double getInverseIntensity(double v) {
        for (int i = 0; i < spline.gridPoints; i++) {
            double rate = spline.I[i];
            if (ratioIsSpline)
                rate *= infectedToRhoSpline.I[i];
            else
                rate *= InfectedToRho.getValue();
            
            v -= rate * (spline.time[i + 1] - spline.time[i]);
            if (v<0){
                v += rate * (spline.time[i + 1] - spline.time[i]);
                // solve for the final time
                return spline.time[i] + v/rate;
            }
        }
		if (ratioIsSpline)
			return spline.time[spline.gridPoints - 1] + v / (spline.I[spline.gridPoints - 1] * infectedToRhoSpline.I[spline.gridPoints - 1]);
		else
	        return spline.time[spline.gridPoints-1] + v/(spline.I[spline.gridPoints-1]*InfectedToRho.getValue());
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

    @Override
    public void init(PrintStream printStream) {
    	if (hasNeSpline)
	        for (int i = 0; i < spline.gridPoints; i+=1) {
	            printStream.print("reassortment" + i + "\t");
	        }
    	else
	        for (int i = 0; i < infectedToRhoSpline.gridPoints; i+=1) {
	            printStream.print("reassortment" + i + "\t");
	        }
    }

    @Override
    public void log(long l, PrintStream printStream) {
    	if (hasNeSpline)
	        for (int i = 0; i < spline.gridPoints; i+=1) {
	        	if (ratioIsSpline)
	        		printStream.print(spline.I[i]*infectedToRhoSpline.I[i] + "\t");
	        	else        		
	        		printStream.print(spline.I[i]*InfectedToRho.getValue() + "\t");
	        }
    	else
            for (int i = 0; i < infectedToRhoSpline.gridPoints; i+=1) {
            		printStream.print(infectedToRhoSpline.I[i] + "\t");
            }

    }

    @Override
    public void close(PrintStream printStream) {

    }


}