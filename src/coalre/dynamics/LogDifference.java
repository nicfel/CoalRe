package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;


@Description("calculates the differences between the entries of a vector")
public class LogDifference extends CalculationNode implements Function {
    final public Input<Function> functionInput = new Input<>("arg", "argument for which the differences for entries is calculated", Validate.REQUIRED);

    enum Mode {integer_mode, double_mode}

    Mode mode;

    boolean needsRecompute = true;
    double expMean;
    double storedExpMean;

    @Override
    public void initAndValidate() {
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return expMean;
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
    	expMean = 0;
	    for (int i = 1; i < functionInput.get().getDimension(); i++) {
	    	expMean += functionInput.get().getArrayValue(i);
	    }
	    expMean /= (functionInput.get().getDimension());       
        expMean = expMean;
        needsRecompute = false;
    }

    @Override
    public double getArrayValue(int dim) {
        if (needsRecompute) {
            compute();
        }
        return expMean;
    }

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
		storedExpMean = expMean;
		super.store();
    }

    @Override
    public void restore() {
		double tmp = storedExpMean;
		storedExpMean = expMean;
		expMean = tmp;
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }
} // class Sum
