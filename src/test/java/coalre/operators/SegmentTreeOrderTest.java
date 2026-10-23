package coalre.operators;

import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import coalre.CoalReTestClass;
import coalre.network.Network;
import coalre.network.NetworkNode;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

public class SegmentTreeOrderTest extends CoalReTestClass {

    /**
     * The network names its segments in the order of its own segmentTree input, and
     * operators map segments to their segment trees by those names. When an operator
     * lists the trees in a different order, the segment trees it rebuilds on its first
     * move must still be the ones its later moves address by name.
     */
    @Test
    public void testOperatorTreeOrderDiffersFromNetwork() {
        Randomizer.setSeed(1);

        TraitSet traitSet = getContempDateTraitSet(getTaxonSet(10));
        List<Tree> trees = getSegmentTreeObjects(3, traitSet);
        Network network = getContempNetwork(trees, 2.0, traitSet);

        List<Tree> operatorTrees = new ArrayList<>(List.of(trees.get(0), trees.get(2), trees.get(1)));

        // Initialise the trees by position in the operator's list, as
        // SegmentTreeInitializer does with the segmentIndex BEAUti assigns.
        for (int j = 0; j < operatorTrees.size(); j++)
            network.updateSegmentTree(operatorTrees.get(j), j);

        UniformNetworkNodeHeightOperator operator = new UniformNetworkNodeHeightOperator();
        operator.initByName("network", network, "segmentTree", operatorTrees, "weight", 1.0);

        // Most operators index segmentTrees directly by segment.
        for (int seg = 0; seg < network.getSegmentCount(); seg++)
            Assert.assertSame(trees.get(seg), operator.segmentTrees.get(seg));

        Assert.assertTrue(network.segmentTreesNeedInit);
        operator.proposal();
        Assert.assertFalse(network.segmentTreesNeedInit);

        for (NetworkNode node : network.getNodes()) {
            for (int seg = 0; seg < network.getSegmentCount(); seg++) {
                if (node.isCoalescence(seg))
                    Assert.assertEquals("segment " + seg, node.getHeight(),
                            operator.getSegmentTreeNode(node, seg).getHeight(), 1e-12);
            }
        }
    }
}
