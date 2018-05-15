package coalre.network;

import beast.core.StateNode;
import coalre.distribution.NetworkEvent;
import coalre.network.parser.NetworkBaseVisitor;
import coalre.network.parser.NetworkLexer;
import coalre.network.parser.NetworkParser;
import coalre.network.parser.NetworkVisitor;
import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTree;
import org.w3c.dom.Node;

import java.io.PrintStream;
import java.util.*;

public class Network extends StateNode {

    protected NetworkEdge rootEdge;
    protected int nSegments;


    public Network() {
    }

    public int getSegmentCount() {
        return nSegments;
    }

    @Override
    public void initAndValidate() { }

	public NetworkEdge getRootEdge() {
        return rootEdge;
    }

    public void setRootEdge(NetworkEdge rootEdge) {
        this.rootEdge = rootEdge;
    }

    public String getExtendedNewick() {
        return rootEdge.getExtendedNewick();
    }

    public Set<NetworkNode> getNodes() {
        Set<NetworkNode> networkNodeSet = new HashSet<>();

        getNodesRecurse(rootEdge, networkNodeSet);

        return networkNodeSet;
    }

    private void getNodesRecurse(NetworkEdge lineage, Set<NetworkNode> networkNodeSet) {

        if (networkNodeSet.contains(lineage.getChildNode()))
            return;

        networkNodeSet.add(lineage.getChildNode());

        for (NetworkEdge childLineage : lineage.getChildNode().getChildEdges())
            getNodesRecurse(childLineage, networkNodeSet);
    }

    public void fromExtendedNewick(String newickStr) {

        CharStream inputStream = CharStreams.fromString(newickStr);
        NetworkLexer lexer = new NetworkLexer(inputStream);
        CommonTokenStream tokenStream = new CommonTokenStream(lexer);
        NetworkParser parser = new NetworkParser(tokenStream);
        ParseTree tree = parser.network();
//        System.out.println(tree.toStringTree(parser));

        NetworkBuilderVisitor builder = new NetworkBuilderVisitor();
        rootEdge = builder.visit(tree);

        System.out.print(toString());
    }

    /** StateNode implementation: **/

    @Override
	public String toString() {
        return getExtendedNewick();
    }

    @Override
    public void setEverythingDirty(boolean isDirty) {
    }

    @Override
    public StateNode copy() {
        return null;
    }

    @Override
    public void assignTo(StateNode other) {

    }

    @Override
    public void assignFrom(StateNode other) {

    }

    @Override
    public void assignFromFragile(StateNode other) {

    }

    @Override
    public void fromXML(Node node) {

    }

    @Override
    public int scale(double scale) {
        return 0;
    }

    @Override
    protected void store() {

    }

    @Override
    public void restore() {

    }

    @Override
    public int getDimension() {
        return 0;
    }

    @Override
    public double getArrayValue() {
        return 0;
    }

    @Override
    public double getArrayValue(int dim) {
        return 0;
    }

    /** Logable implementation: **/

    @Override
    public void init(PrintStream out) {
        out.println("#nexus");
        out.println("begin trees;");
    }

    @Override
    public void close(PrintStream out) {
        out.println("end trees;");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.println("tree STATE_" + sample + " = " + getExtendedNewick());
    }


    /**
     * Visitor class used to build network from parser-generated AST.
     */
    class NetworkBuilderVisitor extends NetworkBaseVisitor<NetworkEdge> {

        Map<Integer, NetworkNode> seenHybrids;
        Map<NetworkEdge, Double> edgeLengths;

        private void convertEdgeLengthsToNodeHeights() {

        }

        double getMaxRootToLeafTime(NetworkNode node, Set<NetworkNode> seenNodes) {

            if (seenNodes.contains(node))
                return 0.0;

            seenNodes.add(node);

            double maxTime = 0.0;
            for (NetworkEdge childEdge : node.getChildEdges()) {
                NetworkNode childNode = childEdge.getChildNode();
                childNode.setHeight(node.getHeight()-edgeLengths.get(childEdge));

                double thisTime = edgeLengths.get(childEdge) +
                        getMaxRootToLeafTime(childNode, seenNodes);
                if (thisTime > maxTime)
                    maxTime = thisTime;
            }

            return maxTime;
        }

        void shiftNodeHeights(double maxRTLT, NetworkNode node, Set<NetworkNode> seenNodes) {
            if (seenNodes.contains(node))
                return;

            seenNodes.add(node);

            node.setHeight(node.getHeight() + maxRTLT);

            for (NetworkEdge childEdge : node.getChildEdges())
                shiftNodeHeights(maxRTLT, childEdge.getChildNode(), seenNodes);
        }

        @Override
        public NetworkEdge visitNetwork(NetworkParser.NetworkContext ctx) {
            seenHybrids = new HashMap<>();
            edgeLengths = new HashMap<>();

            NetworkEdge rootEdge = visit(ctx.node());

            Set<NetworkNode> seenNodes = new HashSet<>();
            NetworkNode rootNode = rootEdge.getChildNode();
            rootNode.setHeight(0.0);
            double maxRTLT = getMaxRootToLeafTime(rootNode, seenNodes);

            seenNodes.clear();
            shiftNodeHeights(maxRTLT, rootEdge.getChildNode(), seenNodes);

            return rootEdge;
        }

        @Override
        public NetworkEdge visitNode(NetworkParser.NodeContext ctx) {

            visit(ctx.post());

            NetworkNode node;

            if (ctx.post().hybrid() != null) {
                int hybridID = Integer.valueOf(ctx.post().hybrid().id.getText());

                if (seenHybrids.containsKey(hybridID)) {
                    node = seenHybrids.get(hybridID);
                } else {
                    node = new NetworkNode();
                    seenHybrids.put(hybridID, node);
                }
            } else {
                node = new NetworkNode();
            }

            if (ctx.post().label() != null)
                node.setTaxonLabel(ctx.post().label().getText());

            for (NetworkParser.NodeContext childNodeCtx : ctx.node()) {
                NetworkEdge childEdge = visit(childNodeCtx);
                childEdge.setParentNode(node);
                node.addChildEdge(childEdge);
            }

            boolean segmentsProcessed = false;
            BitSet hasSegments = new BitSet();
            if (ctx.post().meta() != null
                    || ctx.post().meta().attrib() != null) {

                for (NetworkParser.AttribContext attribCtx : ctx.post().meta().attrib()) {
                    if (!attribCtx.attribKey.getText().equals("segments"))
                        continue;

                    if (attribCtx.attribValue().vector() == null)
                        continue;

                    for (NetworkParser.AttribValueContext attribValueCtx : attribCtx.attribValue().vector().attribValue())
                        hasSegments.set(Integer.valueOf(attribValueCtx.getText()));

                    segmentsProcessed = true;
                    break;
                }

            }

            if (!segmentsProcessed) {
                throw new RuntimeException("Segment attribute missing/malformed " +
                        "for edge in input network string.");
            }

            NetworkEdge edge = new NetworkEdge(null, node, hasSegments);
            node.addParentEdge(edge);

            if (ctx.post().length == null) {
                throw new RuntimeException("Edge missing length in input " +
                        "network string.");
            }

            edgeLengths.put(edge, Double.valueOf(ctx.post().length.getText()));

            return edge;
        }
    }

    /**
     * Main method for debugging.
     *
     * @param args unused
     */
    public static void main(String[] args) {

        Network network = new Network();
        network.fromExtendedNewick("((((t4[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.215550515569861,t3[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.31555051556986097)[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.19210197814368535,(t1[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.5324769040894228,#H0[&segments={2, 5, 6, 7}]:0.18961198047244393)[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.1751755896241235)[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.18917632524163308,t5[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.4968288189551794)[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.19990298145199947,((t10[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.06464178732090486,(t9[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.0384994570159497,t7[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.23849945701594968)[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.1261423303049552)[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.12164810442847562,((t8[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.1428649236169789)#H0[&segments={0, 1, 3, 4}]:0.10219993034575986,(t2[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.06745859531190146,t6[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.16745859531190146)[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.2776062586508373)[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.14122503778664175)[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.5104419086577984)[&segments={0, 1, 2, 3, 4, 5, 6, 7}]:0.0;");
    }
}