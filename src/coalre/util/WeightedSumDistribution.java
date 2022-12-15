package coalre.util;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.util.Randomizer;



@Description("Dirichlet distribution.  p(x_1,...,x_n;alpha_1,...,alpha_n) = 1/B(alpha) prod_{i=1}^K x_i^{alpha_i - 1} " +
        "where B() is the beta function B(alpha) = prod_{i=1}^K Gamma(alpha_i)/ Gamma(sum_{i=1}^K alpha_i}. ")
public class WeightedSumDistribution extends ParametricDistribution {
    public final Input<List<ParametricDistribution>> distInput = new Input<>("distr",
            "distribution used to calculate prior over MRCA time, "
                    + "e.g. normal, beta, gamma. If not specified, monophyletic must be true", new ArrayList<>());
	
    final public Input<RealParameter> weightsInput = new Input<>("weights", "weighting of the individual distrbutions ", Validate.REQUIRED);

    RealParameter weights;
    
    List<ParametricDistribution> distributions;
    
    @Override
    public void initAndValidate() {
    	weights = weightsInput.get();
    	distributions = distInput.get();
    	
    	if (weights.getDimension()!=distributions.size())
    		throw new IllegalArgumentException("the number of weights given differs from the number of distribution");
    	    	
    }

    @Override
    public Distribution getDistribution() {
        return null;
    }

    @Override
    public double calcLogP(Function pX) {
        double logP = 0;
        for (int i = 0; i < pX.getDimension(); i++) {
            double x = pX.getArrayValue(i);
            double prob = 0;
            for (int j = 0; j < weights.getDimension(); j++){
            	try {
					prob += weights.getArrayValue(j)*distributions.get(j).cumulativeProbability(x);
				} catch (MathException e) {
					System.out.println(e);
				}
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
        for (int j = 0; j < weights.getDimension(); j++){
			prob += weights.getArrayValue(j)*distributions.get(j).density(x);
        }
        logP += Math.log(prob);
        return logP;
    }
    
	@Override
	public Double[][] sample(int size) {
		return null;
	}
}
