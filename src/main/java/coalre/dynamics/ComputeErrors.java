package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.Real;
import beast.base.spec.type.RealVector;

import java.util.ArrayList;
import java.util.List;


@Description("calculates the differences between the entries of a vector")
public class ComputeErrors extends CalculationNode implements RealVector<Real> {
    final public Input<RealVector<Real>> functionInput = new Input<>("arg", "argument for which the differences for entries is calculated", Validate.REQUIRED);
    final public Input<RealVector<Real>> casesInput = new Input<>("logCases", "log of the cases", Validate.REQUIRED);
    final public Input<RealVector<Real>> overallNeScalerInput = new Input<>("overallNeScaler", "argument for which the differences for entries is calculated", Validate.REQUIRED);

    enum Mode {integer_mode, double_mode}

    Mode mode;

    boolean needsRecompute = true;
    double[] errorTerm;
    double[] storedErrorTerm;

    @Override
    public void initAndValidate() {
    	errorTerm = new double[functionInput.get().size()];
    	storedErrorTerm = new double[functionInput.get().size()];
    }

    @Override
    public Real getDomain() {
        return Real.INSTANCE;
    }

    @Override
    public int size() {
        return errorTerm.length;
    }

    @Override
    public List<Double> getElements() {
        if (needsRecompute) {
            compute();
        }
        List<Double> elements = new ArrayList<>(errorTerm.length);
        for (double v : errorTerm) elements.add(v);
        return elements;
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {

		for (int i = 0; i < functionInput.get().size(); i++) {
			errorTerm[i] = functionInput.get().get(i) - casesInput.get().get(i) - overallNeScalerInput.get().get(i);
		}
        needsRecompute = false;
    }

    @Override
    public double get(int dim) {
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
