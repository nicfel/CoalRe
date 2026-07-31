package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.Real;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;


@Description("calculates the differences between the entries of a vector")
public class ExpMean extends CalculationNode implements RealScalar<Real> {
    final public Input<RealVector<Real>> functionInput = new Input<>("arg", "argument for which the differences for entries is calculated", Validate.REQUIRED);

    enum Mode {integer_mode, double_mode}

    Mode mode;

    boolean needsRecompute = true;
    double expMean;
    double storedExpMean;

    @Override
    public void initAndValidate() {
    }

    @Override
    public Real getDomain() {
        return Real.INSTANCE;
    }

    @Override
    public double get() {
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
	    for (int i = 1; i < functionInput.get().size(); i++) {
	    	expMean += Math.exp(functionInput.get().get(i));
	    }
	    expMean /= (functionInput.get().size());
	    expMean = Math.exp(expMean);
        needsRecompute = false;
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
