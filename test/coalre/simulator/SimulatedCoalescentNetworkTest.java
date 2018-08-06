package coalre.simulator;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.Randomizer;
import coalre.CoalReTestClass;
import coalre.statistics.NetworkStatsLogger;
import org.junit.Assert;
import org.junit.Test;
import test.beast.beast2vs1.trace.DiscreteStatistics;

import java.util.List;

public class SimulatedCoalescentNetworkTest extends CoalReTestClass {

    @Test
    public void testSimulator() {
        Randomizer.setSeed(1);

        TraitSet dateTrait = getContempDateTraitSet(getTaxonSet(10));
        List<Tree> segmentTrees = getSegmentTreeObjects(8, dateTrait);

        ConstantPopulation populationFunction = new ConstantPopulation();
        populationFunction.initByName("popSize", new RealParameter("1.0"));

        int N = 10000;
        double[] reassortmentNodeCounts = new double[N];
        double[] networkHeights = new double[N];
        double[] networkLengths = new double[N];

        for (int i = 0; i < N; i++) {

            SimulatedCoalescentNetwork network = new SimulatedCoalescentNetwork();
            network.initByName(
                    "reassortmentRate", new RealParameter("1.0"),
                    "populationModel", populationFunction,
                    "segmentTree", segmentTrees.get(0),
                    "segmentTree", segmentTrees.get(1),
                    "segmentTree", segmentTrees.get(2),
                    "segmentTree", segmentTrees.get(3),
                    "segmentTree", segmentTrees.get(4),
                    "segmentTree", segmentTrees.get(5),
                    "segmentTree", segmentTrees.get(6),
                    "segmentTree", segmentTrees.get(7),
                    "enableSegmentTreeUpdate", false);

            reassortmentNodeCounts[i] = NetworkStatsLogger.getReassortmentCount(network);
            networkHeights[i] = NetworkStatsLogger.getTotalHeight(network);
            networkLengths[i] = NetworkStatsLogger.getTotalEdgeLength(network);
        }

        double meanCount = DiscreteStatistics.mean(reassortmentNodeCounts);
        double meanHeight = DiscreteStatistics.mean(networkHeights);
        double meanLength = DiscreteStatistics.mean(networkLengths);

        System.out.println(meanCount);
        System.out.println(meanHeight);
        System.out.println(meanLength);

        Assert.assertEquals(8.05, meanCount, 0.1);
        Assert.assertEquals(2.93, meanHeight, 0.1);
        Assert.assertEquals(10.21, meanLength, 0.5);
    }
}
