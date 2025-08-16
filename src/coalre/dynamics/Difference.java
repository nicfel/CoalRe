package coalre.dynamics;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;


@Description("calculates the differences between the entries of a vector")
public class Difference extends CalculationNode implements Function {
    final public Input<Function> functionInput = new Input<>("arg", "argument for which the differences for entries is calculated", Validate.REQUIRED);
    final public Input<RealParameter> rateShiftInput = new Input<>("rateShift", "rate shift parameter");
    final public Input<Integer> independentAfter = new Input<>("independentAfter", "ignore difference after that index");
    final public Input<Boolean> skipLastInput = new Input<>("skipLast", "skip last element", false);


    boolean needsRecompute = true;
    double[] difference;
    double[] storedDifference;
    int from;
    int to;
    int skipLast = 0;
    

    @Override
    public void initAndValidate() {
		if (skipLastInput.get()) {
			skipLast = 1;
		}
    	if (independentAfter.get()!=null) {
    		difference = new double[functionInput.get().getDimension()-2-skipLast];
    		storedDifference = new double[functionInput.get().getDimension()-2-skipLast];

//    	}else if (independentAfterVector.get()!=null) {
//    		difference = new double[functionInput.get().getDimension()-1 - independentAfterVector.get().getDimension()];
//    		storedDifference = new double[functionInput.get().getDimension()-1 - independentAfterVector.get().getDimension()];
//    		independentAfterList = new ArrayList<>();
//    		for (int i = 0; i < independentAfterVector.get().getDimension(); i++) {
//    			independentAfterList.add((int) independentAfterVector.get().getArrayValue(i));
//    		}
    	}else {
	    	difference = new double[functionInput.get().getDimension()-1];
	    	storedDifference = new double[functionInput.get().getDimension()-1];
    	}
    }

    @Override
    public int getDimension() {
        return difference.length;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return difference[0];
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
    	int offset = 1;
    	if (rateShiftInput.get()==null) {
	        for (int i = 1; i < functionInput.get().getDimension()-skipLast; i++) {
				if (independentAfter.get() != null && i == independentAfter.get()) {
					offset++;
					continue;
				}	        			
	        	difference[i-offset] = functionInput.get().getArrayValue(i-1)-functionInput.get().getArrayValue(i);
	        }    	
    		
    	}else {
	        for (int i = 1; i < functionInput.get().getDimension(); i++) {
				if (independentAfter.get() != null && i == independentAfter.get()) {
					offset++;
					continue;
				} 
	        	difference[i-offset] = (functionInput.get().getArrayValue(i-1)-functionInput.get().getArrayValue(i))
	        			/ (rateShiftInput.get().getArrayValue(i)-rateShiftInput.get().getArrayValue(i-1));
	        }		            		
    	}   	
    	
        
        needsRecompute = false;
    }

    @Override
    public double getArrayValue(int dim) {
        if (needsRecompute) {
            compute();
        }
        return difference[dim];
    }

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
    	System.arraycopy(difference, 0, storedDifference, 0, difference.length);
        super.store();
    }

    @Override
    public void restore() {
    	double [] tmp = storedDifference;
    	storedDifference = difference;
    	difference = tmp;
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }
} // class Sum
