package coalre.util;

import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.distribution.Normal;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;

public class WeightedSumDistributionTest {

    private static Normal normal(double mean, double sigma) {
        Normal normal = new Normal();
        normal.initByName("mean", new RealScalarParam<>(mean, Real.INSTANCE),
                "sigma", new RealScalarParam<>(sigma, PositiveReal.INSTANCE));
        return normal;
    }

    private static Uniform uniform(double lower, double upper) {
        Uniform uniform = new Uniform();
        uniform.initByName("lower", new RealScalarParam<>(lower, Real.INSTANCE),
                "upper", new RealScalarParam<>(upper, Real.INSTANCE));
        return uniform;
    }

    private static WeightedSumDistribution mixture(List<?> components, double... weights) {
        WeightedSumDistribution mixture = new WeightedSumDistribution();
        mixture.initByName("distr", components,
                "weights", new RealVectorParam<>(weights, PositiveReal.INSTANCE));
        return mixture;
    }

    @Test
    public void testDensityAndCdfAreWeightedSums() {
        Normal a = normal(0.0, 1.0);
        Normal b = normal(5.0, 2.0);
        WeightedSumDistribution mixture = mixture(List.of(a, b), 0.3, 0.7);

        for (double x : new double[] {-2.0, 0.0, 1.5, 5.0, 9.0}) {
            Assert.assertEquals(0.3*a.density(x) + 0.7*b.density(x), mixture.density(x), 1e-12);
            Assert.assertEquals(Math.log(0.3*a.density(x) + 0.7*b.density(x)), mixture.logDensity(x), 1e-12);
            Assert.assertEquals(0.3*a.cumulativeProbability(x) + 0.7*b.cumulativeProbability(x),
                    mixture.cumulativeProbability(x), 1e-12);
        }
    }

    @Test
    public void testBoundsSpanAllComponents() {
        Uniform low = uniform(0.0, 1.0);
        Uniform high = uniform(3.0, 4.0);
        WeightedSumDistribution mixture = mixture(List.of(low, high), 1.0, 1.0);

        Assert.assertEquals(0.0, mixture.getLowerBoundOfParameter(), 0.0);
        Assert.assertEquals(4.0, mixture.getUpperBoundOfParameter(), 0.0);
    }

    @Test
    public void testSamplesFollowWeights() {
        Randomizer.setSeed(1);
        Uniform low = uniform(0.0, 1.0);
        Uniform high = uniform(3.0, 4.0);
        WeightedSumDistribution mixture = mixture(List.of(low, high), 1.0, 3.0);

        int n = 20000, nHigh = 0;
        for (int i = 0; i < n; i++) {
            double x = mixture.sample().get(0);
            Assert.assertTrue((x >= 0.0 && x <= 1.0) || (x >= 3.0 && x <= 4.0));
            if (x >= 3.0)
                nHigh++;
        }
        Assert.assertEquals(0.75, nHigh/(double) n, 0.02);
    }
}
