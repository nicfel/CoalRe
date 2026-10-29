package coalre.util;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;
import beast.base.util.Randomizer;



@Description("Weighted sum of scalar distributions: density and CDF are the weighted sums of "
        + "those of the component distributions.")
public class WeightedSumDistribution extends ScalarDistribution<RealScalar<Real>, Double> {
    public final Input<List<ScalarDistribution>> distInput = new Input<>("distr",
            "distribution used to calculate prior over MRCA time, "
                    + "e.g. normal, beta, gamma. If not specified, monophyletic must be true", new ArrayList<>());

    final public Input<RealVector<PositiveReal>> weightsInput = new Input<>("weights", "weighting of the individual distrbutions ", Validate.REQUIRED);

    RealVector<PositiveReal> weights;

    List<ScalarDistribution> distributions;

    @Override
    public void initAndValidate() {
    	weights = weightsInput.get();
    	distributions = distInput.get();

    	if (weights.size()!=distributions.size())
    		throw new IllegalArgumentException("the number of weights given differs from the number of distribution");

        super.initAndValidate();
    }

    @Override
    public void refresh() {
        // nothing to cache: each component refreshes itself when evaluated
    }

    // Scores the parameter by the weighted sum of the component CDFs, not their densities.
    // This is the BEAST 2 behaviour (calcLogP(Function)), kept so that the model is unchanged.
    @Override
    protected double calcLogP(Double x) {
        return Math.log(cumulativeProbability(x));
    }

    @Override
    public double density(double x) {
        double prob = 0;
        for (int j = 0; j < weights.size(); j++){
			prob += weights.get(j)*distributions.get(j).density(x);
        }
        return prob;
    }

    @Override
    public double cumulativeProbability(double x) {
        double prob = 0;
        for (int j = 0; j < weights.size(); j++){
            prob += weights.get(j)*distributions.get(j).cumulativeProbability(x);
        }
        return prob;
    }

    @Override
    public List<Double> sample() {
        double totalWeight = 0;
        for (int j = 0; j < weights.size(); j++)
            totalWeight += weights.get(j);

        double u = Randomizer.nextDouble()*totalWeight;
        int j = 0;
        while (j < weights.size()-1 && u >= weights.get(j)) {
            u -= weights.get(j);
            j++;
        }
        return List.of(((Number) distributions.get(j).sample().get(0)).doubleValue());
    }

    @Override
    public Double getLowerBoundOfParameter() {
        double lower = Double.POSITIVE_INFINITY;
        for (ScalarDistribution dist : distributions)
            lower = Math.min(lower, ((Number) dist.getLowerBoundOfParameter()).doubleValue());
        return lower;
    }

    @Override
    public Double getUpperBoundOfParameter() {
        double upper = Double.NEGATIVE_INFINITY;
        for (ScalarDistribution dist : distributions)
            upper = Math.max(upper, ((Number) dist.getUpperBoundOfParameter()).doubleValue());
        return upper;
    }
}
