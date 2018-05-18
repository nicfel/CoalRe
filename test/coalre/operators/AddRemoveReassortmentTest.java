package coalre.operators;

import coalre.CoalReTestClass;
import coalre.network.Network;
import coalre.network.NetworkNode;
import org.junit.Test;

import java.util.*;

public class AddRemoveReassortmentTest extends CoalReTestClass {

    String networkString = "((#H0[&\"segments\"={1,2,4,7}]:0.09623670174825327," +
            "\"t5\"[&\"segments\"={0,1,2,3,4,5,6,7}]:0.9532676725160346)" +
            "[&\"segments\"={0,1,2,3,4,5,6,7}]:0.43959702603952655,(((#H1" +
            "[&\"segments\"={0,7}]:0.5444373856629039,((\"t3\"" +
            "[&\"segments\"={0,1,2,3,4,5,6,7}]:0.27262918122452245,#H2" +
            "[&\"segments\"={1,3,6}]:0.05445266311806374)" +
            "[&\"segments\"={0,1,2,3,4,5,6,7}]:0.12727797841969501,\"t4\"" +
            "[&\"segments\"={0,1,2,3,4,5,6,7}]:0.2999071596442175)" +
            "[&\"segments\"={0,1,2,3,4,5,6,7}]:0.4250628592896514)" +
            "[&\"segments\"={0,1,2,3,4,5,6,7}]:0.09731951917853188,(((\"t1\"" +
            "[&\"segments\"={0,1,2,3,4,5,6,7}]:0.41817651810645873)#H2" +
            "[&\"segments\"={0,2,4,5,7}]:0.06235611516450629)#H1" +
            "[&\"segments\"={2,4,5}]:0.10585974103145485,\"t2\"" +
            "[&\"segments\"={0,1,2,3,4,5,6,7}]:0.4863923743024199)" +
            "[&\"segments\"={0,1,2,3,4,5,6,7}]:0.5358971638099809)" +
            "[&\"segments\"={0,1,2,3,4,5,6,7}]:0.13474143265538063)#H0" +
            "[&\"segments\"={0,3,5,6}]:0.5358337277877798)" +
            "[&\"segments\"={0,1,2,3,4,5,6,7}]:0.0;";

    @Test
    public void testAddSegment() {

        Network network = new Network(networkString);

        NetworkNode leafNode = new ArrayList<>(network.getLeafNodes()).get(0);

        AddRemoveReassortment operator = new AddRemoveReassortment();

        BitSet segmentsToAdd = new BitSet();
        segmentsToAdd.set(50);
        operator.addSegmentsToAncestors(leafNode.getParentEdges().get(0),
                segmentsToAdd, new HashSet<>());

        System.out.println(network);
    }

    @Test
    public void testRemoveSegment() {

        Network network = new Network(networkString);

        NetworkNode leafNode = new ArrayList<>(network.getLeafNodes()).get(0);

        AddRemoveReassortment operator = new AddRemoveReassortment();

        BitSet segmentsToRemove = new BitSet();
        segmentsToRemove.set(0);
        operator.removeSegmentsFromAncestors(leafNode.getParentEdges().get(0),
                segmentsToRemove, new HashSet<>());

        System.out.println(network);
    }
}
