package coalre.util;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.type.RealVector;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.util.Randomizer;



@Description("Dirichlet distribution.  p(x_1,...,x_n;alpha_1,...,alpha_n) = 1/B(alpha) prod_{i=1}^K x_i^{alpha_i - 1} " +
        "where B() is the beta function B(alpha) = prod_{i=1}^K Gamma(alpha_i)/ Gamma(sum_{i=1}^K alpha_i}. ")
public class WeightedSumDistribution extends ParametricDistribution {
    public final Input<List<ParametricDistribution>> distInput = new Input<>("distr",
            "distribution used to calculate prior over MRCA time, "
                    + "e.g. normal, beta, gamma. If not specified, monophyletic must be true", new ArrayList<>());
	
    final public Input<RealVector<PositiveReal>> weightsInput = new Input<>("weights", "weighting of the individual distrbutions ", Validate.REQUIRED);

    RealVector<PositiveReal> weights;
    
    List<ParametricDistribution> distributions;
    
    @Override
    public void initAndValidate() {
    	weights = weightsInput.get();
    	distributions = distInput.get();
    	
    	if (weights.size()!=distributions.size())
    		throw new IllegalArgumentException("the number of weights given differs from the number of distribution");
    	    	
    }

    // BEAST3: ParametricDistribution.getDistribution() returns Object
    // (commons-math2 Distribution was replaced by commons-statistics).
    @Override
    public Object getDistribution() {
        return null;
    }

    // BEAST3: override the spec RealVector overload — with typed inputs the
    // legacy calcLogP(Function) overload is never called.
    @Override
    public double calcLogP(RealVector<?> pX) {
        double logP = 0;
        for (int i = 0; i < pX.size(); i++) {
            double x = pX.get(i);
            double prob = 0;
            for (int j = 0; j < weights.size(); j++){
                // BEAST3: cumulativeProbability no longer throws MathException.
                prob += weights.get(j)*distributions.get(j).cumulativeProbability(x);
            }
            logP += Math.log(prob);
        }
        return logP;
    }

    @Override
    public double logDensity(double val){
        double logP = 0;
        double x = val;
        double prob = 0;
        for (int j = 0; j < weights.size(); j++){
			prob += weights.get(j)*distributions.get(j).density(x);
        }
        logP += Math.log(prob);
        return logP;
    }
    
	@Override
	public Double[][] sample(int size) {
		return null;
	}
}
