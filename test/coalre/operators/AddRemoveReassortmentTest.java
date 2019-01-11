package coalre.operators;

import beast.util.Randomizer;
import coalre.CoalReTestClass;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import org.junit.Assert;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

public class AddRemoveReassortmentTest extends CoalReTestClass {

    String networkString = "((#H0[&segments={1, 2, 4, 7},segsCarried=4]:0.09623670174825327,t5" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.9532676725160347)" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.43959702603952655,(((#H1" +
            "[&segments={0, 7},segsCarried=2]:0.5444373856629039,((t3" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.2726291812245225,#H2" +
            "[&segments={1, 3, 6},segsCarried=3]:0.05445266311806374)" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.1272779784196949,t4" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.29990715964421755)" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.4250628592896515)" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.09731951917853188,(((t1" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.41817651810645873)#H2" +
            "[&segments={0, 2, 4, 5, 7},segsCarried=5]:0.06235611516450623)#H1" +
            "[&segments={2, 4, 5},segsCarried=3]:0.10585974103145479,t2" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.4863923743024199)" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.535897163809981)" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.13474143265538063)#H0" +
            "[&segments={0, 3, 5, 6},segsCarried=4]:0.5358337277877798)" +
            "[&segments={0, 1, 2, 3, 4, 5, 6, 7},segsCarried=8]:0.0;";

    @Test
    public void testAddRemoveSegment() {

        Network network = new Network(networkString);
        NetworkNode leafNode = new ArrayList<>(network.getLeafNodes()).get(0);

        AddRemoveReassortment operator = new AddRemoveReassortment();

        BitSet segmentsToAdd = new BitSet();
        segmentsToAdd.set(8, 15);

        double logPadd = operator.addSegmentsToAncestors(
                leafNode.getParentEdges().get(0), segmentsToAdd);

        BitSet allSegments = new BitSet();
        allSegments.set(0, 15);
        Assert.assertEquals(allSegments, network.getRootEdge().hasSegments);

        BitSet segmentsToRemove = new BitSet();
        segmentsToRemove.set(8, 15);

        double logPremove = operator.removeSegmentsFromAncestors(
                leafNode.getParentEdges().get(0), segmentsToRemove);

        Assert.assertEquals(networkString, network.toString());
        Assert.assertEquals(logPadd, logPremove, 1e-10);
    }

    @Test
    public void testAddRemoveReassortmentEdge() {
        Network network = getContempNetwork(2, 8, 0.0);

        AddRemoveReassortment operator = new AddRemoveReassortment();
        operator.initByName("alpha", 1.0,
                "network", network,
                "weight", 1.0);

        NetworkNode origRoot = network.getRootEdge().childNode;

        NetworkEdge sourceEdge = network.getRootEdge().childNode.getChildEdges().get(0);
        double sourceTime = sourceEdge.getLength()/2.0;
        NetworkEdge destEdge = network.getRootEdge();
        double destTime = destEdge.childNode.getHeight() + 1.0;

        double logP1 = operator.addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);

        NetworkEdge edgeToRemove = sourceEdge.parentNode.getParentEdges().get(0).parentNode == origRoot
                ? sourceEdge.parentNode.getParentEdges().get(1)
                : sourceEdge.parentNode.getParentEdges().get(0);

        double logP2 = operator.removeReassortmentEdge(edgeToRemove);

        Assert.assertEquals(logP1, -logP2, 1e-10);

        sourceEdge = network.getRootEdge().childNode.getChildEdges().get(1);
        sourceTime = sourceEdge.getLength()/4.0;
        destEdge = sourceEdge;
        destTime = sourceEdge.getLength()*3.0/4.0;

        logP1 = operator.addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);

        edgeToRemove = sourceEdge.parentNode.getParentEdges().get(0);

        logP2 = operator.removeReassortmentEdge(edgeToRemove);

        Assert.assertEquals(logP1, -logP2, 1e-10);
    }

    @Test
    public void testRemoveReassortment() {
        // TODO Flesh out this test

        Network network = new Network(networkString);

        AddRemoveReassortment operator = new AddRemoveReassortment();
        operator.initByName("network", network,
                "alpha", 1.0,
                "weight", 1.0);

        System.out.println(network.getExtendedNewickVerbose());

        operator.removeReassortment();

        System.out.println(network.getExtendedNewickVerbose());
    }

    @Test
    public void testAddReassortment() {
        // TODO Flesh out this test

        Network network = new Network(networkString);

        AddRemoveReassortment operator = new AddRemoveReassortment();
        operator.initByName("network", network,
                "alpha", 1.0,
                "weight", 1.0);

        System.out.println(network.getExtendedNewickVerbose());

        double logHR = operator.addReassortment();

        System.out.println(network.getExtendedNewickVerbose());

        System.out.println(logHR);
    }
}
