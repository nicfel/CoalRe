package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;


@Description("calculates the differences between the entries of a vector")
public class ComputeErrors extends CalculationNode implements Function {
    final public Input<Function> functionInput = new Input<>("arg", "argument for which the differences for entries is calculated", Validate.REQUIRED);
    final public Input<Function> casesInput = new Input<>("logCases", "log of the cases", Validate.REQUIRED);
    final public Input<Function> overallNeScalerInput = new Input<>("overallNeScaler", "argument for which the differences for entries is calculated", Validate.REQUIRED);

    enum Mode {integer_mode, double_mode}

    Mode mode;

    boolean needsRecompute = true;
    double[] errorTerm;
    double[] storedErrorTerm;

    @Override
    public void initAndValidate() {
    	errorTerm = new double[functionInput.get().getDimension()];
    	storedErrorTerm = new double[functionInput.get().getDimension()];
    }

    @Override
    public int getDimension() {
        return errorTerm.length;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return errorTerm[0];
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
    	
		for (int i = 0; i < functionInput.get().getDimension(); i++) {
			errorTerm[i] = functionInput.get().getArrayValue(i) - casesInput.get().getArrayValue(i) - overallNeScalerInput.get().getArrayValue(i);
		}
        needsRecompute = false;
    }

    @Override
    public double getArrayValue(int dim) {
        if (needsRecompute) {
            compute();
        }
        return errorTerm[dim];
    }

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
    	System.arraycopy(errorTerm, 0, storedErrorTerm, 0, errorTerm.length);
        super.store();
    }

    @Override
    public void restore() {
    	double [] tmp = storedErrorTerm;
    	storedErrorTerm = errorTerm;
    	errorTerm = tmp;
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }
} // class Sum
