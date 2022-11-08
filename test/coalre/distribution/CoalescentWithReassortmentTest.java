package coalre.distribution;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import coalre.CoalReTestClass;
import coalre.network.Network;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

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
        populationFunction.initByName("popSize", new RealParameter("1.0"));

        CoalescentWithReassortment coalWR = new CoalescentWithReassortment();
        coalWR.initByName("networkIntervals", networkIntervals,
                "reassortmentRate", new RealParameter("1.0"),
                "populationModel", populationFunction);

        assertEquals(-16.258280263919616, coalWR.calculateLogP(), 1e-10);
    }
}
