package coalre.simulator;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.ConstantPopulation;
import coalre.CoalReTestClass;
import org.junit.Test;

import java.util.List;

public class SimulateCoalescentWithReassortmentTest extends CoalReTestClass {

    @Test
    public void testSimulator() {

        /*
        This isn't actually a test yet: just some code used in debugging simulator.
        But this will form the basis for a statistical test of the simulation
        code and a separate test of the serialization code.
         */

        TraitSet dateTrait = getContempDateTraitSet(getTaxonSet(10));
        List<Tree> segmentTrees = getSegmentTreeObjects(8, dateTrait);

        ConstantPopulation populationFunction = new ConstantPopulation();
        populationFunction.initByName("popSize", new RealParameter("1.0"));


        SimulateCoalescentWithReassortment network = new SimulateCoalescentWithReassortment();
        network.initByName(
                "rRate", 1.0,
                "populationModel", populationFunction,
                "segmentTree", segmentTrees.get(0),
                "segmentTree", segmentTrees.get(1),
                "segmentTree", segmentTrees.get(2),
                "segmentTree", segmentTrees.get(3),
                "segmentTree", segmentTrees.get(4),
                "segmentTree", segmentTrees.get(5),
                "segmentTree", segmentTrees.get(6),
                "segmentTree", segmentTrees.get(7));

        System.out.println(network);
    }
}
