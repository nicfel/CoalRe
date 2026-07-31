package coalre.distribution;

import beast.base.spec.evolution.tree.coalescent.ConstantPopulation;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import coalre.CoalReTestClass;
import coalre.network.Network;

import org.junit.Assert;
import org.junit.Test;

public class CoalescentWithReassortmentTest extends CoalReTestClass {

    @Test
    public void testDensity() {

        Network network = new Network(
                "(((A[&segments={0,1,2,3,4,5,6,7}]:1)#1[&segments={1,4}]:1," +
                        "B[&segments={0,1,2,3,4,5,6,7}]:2)[&segments={0,1,2,3,4,5,6,7}]:1," +
                        "#1[&segments={0,2,3,5,6,7}]:2)[&segments={0,1,2,3,4,5,6,7}]:0.0;");

        NetworkIntervals networkIntervals = new NetworkIntervals();
        networkIntervals.initByName("network", network);

        ConstantPopulation populationFunction = new ConstantPopulation();
        populationFunction.initByName("popSize", new RealScalarParam<>(1.0, PositiveReal.INSTANCE));

        CoalescentWithReassortment coalWR = new CoalescentWithReassortment();
        coalWR.initByName("networkIntervals", networkIntervals,
                // BEAST3: reassortmentRate is now a typed spec param, not a Function
                "reassortmentRate", new RealScalarParam<>(1.0, PositiveReal.INSTANCE),
                "populationModel", populationFunction);

        Assert.assertEquals(-16.258280263919616, coalWR.calculateLogP(), 1e-10);
    }
}
