package coalre.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;


@Description("joint operator to keep the total reassortment the same")
public class ChangePredictorOperator extends Operator {
	final public Input<List<RealParameter>> predictorInput = new Input<>("predictor", "predictor parameters that are used to calculate the Ne", new ArrayList<>());
	final public Input<RealParameter> NeToReassortmentInput = new Input<>("neToReassortment",
			"the value that maps the number of infected or the Ne to the reassortment rate ");
	final public Input<IntegerParameter> predictorIsActiveInput = new Input<>("predictorIsActive",
			"indicates which predictors are active at which time point");
	final public Input<Integer> independentAfterInput = new Input<>("independentAfter",
			"ignore differences after that index");
	final public Input<RealParameter> effectSizeInput = new Input<>("effectSize",
			"the effect size of the predictors on the reassortment rates", Input.Validate.REQUIRED);
	
	RealParameter NeToReassortment;
	RealParameter effectSize;

	List<RealParameter> predictors;

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

        IntegerParameter param = (IntegerParameter) InputUtil.get(predictorIsActiveInput, this);

        int oldValue = param.getValue(0);
        int i = Randomizer.nextInt(param.getDimension());
        int newValue = Randomizer.nextInt(param.getUpper() - param.getLower() + 1) + param.getLower();

        param.setValue(i, newValue);
        
        double[] currentRates = calculateRates(oldValue);
        
        // set the new values of NeToReassortment, such that the total reassortment rate remains the same
        double[] newRates = calculateRates(newValue);
        for (int j = 0; j < NeToReassortment.getDimension()-1; j++) {
			NeToReassortment.setValue(j, NeToReassortment.getArrayValue(j) + currentRates[j] - newRates[j]);
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
				mean += predictors.get(predictorIndex).getArrayValue(i);
			}
			mean /= (independentAfterInput.get()+1);
			for (int i = 0; i < independentAfterInput.get() + 1; i++) {
				logStandardPredictor[i] = predictors.get(predictorIndex).getArrayValue(i) - mean;
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
		
		
		double[] rates = new double[NeToReassortment.getDimension()];
		if (predictorIndex<predictors.size())  {			
			for (int i = 0; i < independentAfterInput.get()+1; i++) {
				rates[i] = effectSize.getArrayValue(predictorIndex)*
						logStandardPredictor[i] + NeToReassortment.getArrayValue(i);
			}
			for (int i = independentAfterInput.get()+1; i < NeToReassortment.getDimension(); i++) {
				rates[i] = NeToReassortment.getArrayValue(i);
			}
		}else {
			for (int i = 0; i < NeToReassortment.getDimension(); i++) {
				rates[i] = NeToReassortment.getArrayValue(i);
			}
		}
		return rates;
	}

    @Override
    public void optimize(double logAlpha) {
        // nothing to optimise
    }

} // class IntUniformOperator