package coalre.operators;

import java.text.DecimalFormat;
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
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.IntVector;
import beast.base.spec.type.RealVector;
import beast.base.util.Randomizer;
import cern.colt.Arrays;
import coalre.network.Network;


@Description("joint operator to keep the total reassortment the same")
public class EffectSizePredictorOperator extends Operator {
	final public Input<List<RealVector<Real>>> predictorInput = new Input<>("predictor", "predictor parameters that are used to calculate the Ne", new ArrayList<>());
	final public Input<RealVectorParam<Real>> NeToReassortmentInput = new Input<>("neToReassortment",
			"the value that maps the number of infected or the Ne to the reassortment rate ");
	final public Input<IntVector<NonNegativeInt>> predictorIsActiveInput = new Input<>("predictorIsActive",
			"indicates which predictors are active at which time point");
	final public Input<Integer> independentAfterInput = new Input<>("independentAfter",
			"ignore differences after that index");
	final public Input<RealVectorParam<Real>> effectSizeInput = new Input<>("effectSize",
			"the effect size of the predictors on the reassortment rates", Input.Validate.REQUIRED);

    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<Boolean> gaussianInput = new Input<>("gaussian", "Gaussian (=true=default) or uniform delta", true);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, " +
            "i.e. -1. (when positive, gets multiplied by network-height/log2(n-taxa).", -1.0);


	RealVectorParam<Real> NeToReassortment;
	RealVectorParam<Real> effectSize;

	List<RealVector<Real>> predictors;
	
    // shadows size
    double size;
    private double limit;


    @Override
	public void initAndValidate() {
        NeToReassortment = NeToReassortmentInput.get();
        predictors = predictorInput.get();
        effectSize = effectSizeInput.get();       

        size = sizeInput.get();
        limit = limitInput.get();

    }

    /**
     * override this for proposals,
     * returns log of hastingRatio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
    	
    	
    	int activePredictorIndex = predictorIsActiveInput.get().get(0);
    	if (activePredictorIndex >= predictors.size()) {
    		return Double.NEGATIVE_INFINITY; // no active predictor, no proposal
    	}

        double[] currentRates = calculateRates(activePredictorIndex);
//        System.out.println(Arrays.toString(currentRates));


        // propose a new effect size for activePredictorIndex
        double delta = getDelta();
        effectSize.set(activePredictorIndex, effectSize.get(activePredictorIndex) + delta);


        // set the new values of NeToReassortment, such that the total reassortment rate remains the same
        double[] newRates = calculateRates(activePredictorIndex);
        for (int j = 0; j < NeToReassortment.size()-1; j++) {
			NeToReassortment.set(j, NeToReassortment.get(j) + currentRates[j] - newRates[j]);
		}
//        System.out.println(Arrays.toString(calculateRates(activePredictorIndex)));

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
	
    private double getDelta() {
        if (!gaussianInput.get()) {
            return (Randomizer.nextDouble() * size) - (size / 2.0);
        } else {
            return Randomizer.nextGaussian() * size;
        }
    }


	/**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(size);
            final double f = Math.exp(delta);
            if( limit > 0 ) {
                final double lim = limit;
                if( f <= lim ) {
                    size = f;
                }
            } else {
               size = f;
            }
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return size;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        size = value;
    }
    
    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;

        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        final double newDelta = size * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try decreasing size to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try increasing size to about " + formatter.format(newDelta);
        } else return "";
    }
} // class IntUniformOperator