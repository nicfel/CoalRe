package coalre.network;

import org.junit.Assert;
import org.junit.Test;

import java.util.BitSet;

public class NetworkTest {

    @Test
    public void parserTest() {

        BitSet hasSegments = new BitSet();
        hasSegments.set(0, 8);

        NetworkNode nodeA = new NetworkNode();
        nodeA.setTaxonLabel("A");
        nodeA.setHeight(0.0);
        NetworkNode nodeAprime = new NetworkNode();
        nodeAprime.setHeight(1.0);
        NetworkNode nodeB = new NetworkNode();
        nodeB.setTaxonLabel("B");
        nodeB.setHeight(0.0);
        NetworkNode nodeAB = new NetworkNode();
        nodeAB.setHeight(2.0);
        NetworkNode rootNode = new NetworkNode();
        rootNode.setHeight(3.0);

        NetworkEdge edgeA = new NetworkEdge();
        edgeA.hasSegments = (BitSet)hasSegments.clone();
        nodeA.addParentEdge(edgeA);

        nodeAprime.addChildEdge(edgeA);

        NetworkEdge edgeAprime1 = new NetworkEdge();
        edgeAprime1.hasSegments = new BitSet();
        edgeAprime1.hasSegments.set(1);
        edgeAprime1.hasSegments.set(4);
        nodeAprime.addParentEdge(edgeAprime1);

        NetworkEdge edgeAprime2 = new NetworkEdge();
        edgeAprime2.hasSegments = (BitSet)hasSegments.clone();
        edgeAprime2.hasSegments.xor(edgeAprime1.hasSegments);
        nodeAprime.addParentEdge(edgeAprime2);

        NetworkEdge edgeB = new NetworkEdge();
        edgeB.hasSegments = (BitSet)hasSegments.clone();
        nodeB.addParentEdge(edgeB);

        nodeAB.addChildEdge(edgeAprime1).addChildEdge(edgeB);

        NetworkEdge edgeAB =  new NetworkEdge();
        edgeAB.hasSegments = (BitSet)hasSegments.clone();
        nodeAB.addParentEdge(edgeAB);

        rootNode.addChildEdge(edgeAB);
        rootNode.addChildEdge(edgeAprime2);

        NetworkEdge rootEdge = new NetworkEdge(null, rootNode, hasSegments);
        Network network = new Network(rootEdge);

        Network parsedNetwork = new Network("(((A[&segments={0,1,2,3,4,5,6,7}]:1)#1[&segments={1,4}]:1,B[&segments={0,1,2,3,4,5,6,7}]:2)[&segments={0,1,2,3,4,5,6,7}]:1,#1[&segments={0,2,3,5,6,7}]:2)[&segments={0,1,2,3,4,5,6,7}]:0.0;");

        Assert.assertEquals(network.toString(), parsedNetwork.toString());

        // Note to future self:
        // This test may fail even if the parser is correct, as it assumes the
        // generated newick string is unique: it is not.
        // Root of the problem: testing for network equivalence is hard, as
        // discussed at https://en.wikipedia.org/wiki/Graph_isomorphism_problem
    }
}
