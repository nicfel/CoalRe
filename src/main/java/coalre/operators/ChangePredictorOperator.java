package coalre.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.util.InputUtil;
import beast.base.spec.domain.NonNegativeInt;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealVector;
import beast.base.util.Randomizer;


@Description("joint operator to keep the total reassortment the same")
public class ChangePredictorOperator extends Operator {
	final public Input<List<RealVector<Real>>> predictorInput = new Input<>("predictor", "predictor parameters that are used to calculate the Ne", new ArrayList<>());
	final public Input<RealVectorParam<Real>> NeToReassortmentInput = new Input<>("neToReassortment",
			"the value that maps the number of infected or the Ne to the reassortment rate ");
	final public Input<IntScalarParam<NonNegativeInt>> predictorIsActiveInput = new Input<>("predictorIsActive",
			"index of the active predictor, or the number of predictors for none");
	final public Input<Integer> independentAfterInput = new Input<>("independentAfter",
			"ignore differences after that index");
	final public Input<RealVector<Real>> effectSizeInput = new Input<>("effectSize",
			"the effect size of the predictors on the reassortment rates", Input.Validate.REQUIRED);

	RealVectorParam<Real> NeToReassortment;
	RealVector<Real> effectSize;

	List<RealVector<Real>> predictors;

    @Override
	public void initAndValidate() {
        NeToReassortment = NeToReassortmentInput.get();
        predictors = predictorInput.get();
        effectSize = effectSizeInput.get();

    }

    /**
     * override this for proposals,
     * returns log of hastingRatio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        IntScalarParam<NonNegativeInt> param = (IntScalarParam<NonNegativeInt>) InputUtil.get(predictorIsActiveInput, this);

        int oldValue = param.get();
        // Valid states are the predictor indices 0..predictors.size()-1 plus predictors.size(),
        // which means "no predictor active" (see the guards in GLMReassortmentRates).
        //
        // BEAST 2 read this range off the parameter itself, which carried lower/upper set in the
        // XML (upper="4" for 4 predictors). BEAST 3 derives bounds from the Domain instead, and
        // NonNegativeInt inherits getUpper() == Integer.MAX_VALUE from Int -- so the old
        // expression became nextInt(Integer.MAX_VALUE - 0 + 1), which overflows to
        // Integer.MIN_VALUE. The range is therefore taken from the predictor list directly.
        int newValue = Randomizer.nextInt(predictors.size() + 1);

        param.set(newValue);

        double[] currentRates = calculateRates(oldValue);

        // set the new values of NeToReassortment, such that the total reassortment rate remains the same
        double[] newRates = calculateRates(newValue);
        for (int j = 0; j < NeToReassortment.size()-1; j++) {
			NeToReassortment.set(j, NeToReassortment.get(j) + currentRates[j] - newRates[j]);
		}
        return 0.0;
    }
    
	// computes the Ne's at the break points
	private double[] calculateRates(int predictorIndex) {
		double[] logStandardPredictor = new double[independentAfterInput.get()+1];
		if (predictorIndex<predictors.size())  {	
			double mean = 0.0;
			for (int i = 0; i < predictors.size(); i++) {
			}
			for (int i = 0; i < independentAfterInput.get()+1; i++) {
				mean += predictors.get(predictorIndex).get(i);
			}
			mean /= (independentAfterInput.get()+1);
			for (int i = 0; i < independentAfterInput.get() + 1; i++) {
				logStandardPredictor[i] = predictors.get(predictorIndex).get(i) - mean;
			}
			double sd = 0.0;
			for (int i = 0; i < independentAfterInput.get() + 1; i++) {
				sd += Math.pow(logStandardPredictor[i], 2);
			}
			sd = Math.sqrt(sd / (independentAfterInput.get() + 1));
			for (int i = 0; i < independentAfterInput.get() + 1; i++) {
				logStandardPredictor[i] /= sd;
			}
		}
		
		
		double[] rates = new double[NeToReassortment.size()];
		if (predictorIndex<predictors.size())  {
			for (int i = 0; i < independentAfterInput.get()+1; i++) {
				rates[i] = effectSize.get(predictorIndex)*
						logStandardPredictor[i] + NeToReassortment.get(i);
			}
			for (int i = independentAfterInput.get()+1; i < NeToReassortment.size(); i++) {
				rates[i] = NeToReassortment.get(i);
			}
		}else {
			for (int i = 0; i < NeToReassortment.size(); i++) {
				rates[i] = NeToReassortment.get(i);
			}
		}
		return rates;
	}

    @Override
    public void optimize(double logAlpha) {
        // nothing to optimise
    }

} // class IntUniformOperator