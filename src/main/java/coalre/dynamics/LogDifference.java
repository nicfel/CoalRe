package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.type.RealVector;

import java.util.ArrayList;
import java.util.List;


@Description("calculates the differences between the entries of a vector")
public class LogDifference extends CalculationNode implements RealVector<Real> {
    final public Input<RealVector<Real>> functionInput = new Input<>("arg", "argument for which the differences for entries is calculated", Validate.REQUIRED);
    final public Input<RealVector<NonNegativeReal>> rateShiftInput = new Input<>("rateShift", "rate shift parameter");

    final public Input<Integer> fromInput = new Input<>("from", "start point for diff calculation", 0);
    final public Input<Integer> toInput = new Input<>("to" , "end point for diff calculation",-1);


    boolean needsRecompute = true;
    double[] difference;
    double[] storedDifference;
    int from;
    int to;


    @Override
    public void initAndValidate() {
    	from = fromInput.get();
    	to = toInput.get();
		if (to < 0) {
			to = functionInput.get().size() - 1;
		}


    	difference = new double[to-from+1];
    	storedDifference = new double[to-from+1];
    }

    @Override
    public Real getDomain() {
        return Real.INSTANCE;
    }

    @Override
    public int size() {
        return difference.length;
    }

    @Override
    public List<Double> getElements() {
        if (needsRecompute) {
            compute();
        }
        List<Double> elements = new ArrayList<>(difference.length);
        for (double v : difference) elements.add(v);
        return elements;
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
    	double mean = 0;
        for (int i = from; i <= to; i++) {
        	difference[i-from] = functionInput.get().get(i);
        	mean += difference[i-from];
        }
        mean /= (double) difference.length;
		for (int i = 0; i < difference.length; i++) {
			difference[i] -= mean;
			difference[i] = Math.abs(difference[i])+1e-16; // add small value to avoid log(0);
		}
		needsRecompute = false;

    }

    @Override
    public double get(int dim) {
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
