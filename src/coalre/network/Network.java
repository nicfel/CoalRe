package coalre.network;

import beast.core.StateNode;
import beast.evolution.tree.Tree;
import coalre.network.parser.NetworkBaseVisitor;
import coalre.network.parser.NetworkLexer;
import coalre.network.parser.NetworkParser;
import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTree;
import org.w3c.dom.Node;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

public class Network extends StateNode {

    protected NetworkEdge rootEdge;

    protected NetworkEdge storedRootEdge;

    protected Integer segmentCount = null;

    public Network() {
    }

    public Network(NetworkEdge rootEdge) {
        this.rootEdge = rootEdge;
    }

    public Network(String newickString) {
        fromExtendedNewick(newickString);
    }

    @Override
    public void initAndValidate() { }

    /**
     * @return the root edge of the network
     */
	public NetworkEdge getRootEdge() {
        return rootEdge;
    }

    /**
     * Set the root edge of the network.
     *
     * @param rootEdge the new root edge.
     */
    public void setRootEdge(NetworkEdge rootEdge) {
        this.rootEdge = rootEdge;
    }

    /**
     * @return set of node objects comprising network
     */
    public Set<NetworkNode> getNodes() {
        Set<NetworkNode> networkNodeSet = new HashSet<>();

        getNodesRecurse(rootEdge, networkNodeSet);

        return networkNodeSet;
    }

    private void getNodesRecurse(NetworkEdge lineage, Set<NetworkNode> networkNodeSet) {

        if (networkNodeSet.contains(lineage.childNode))
            return;

        networkNodeSet.add(lineage.childNode);

        for (NetworkEdge childLineage : lineage.childNode.getChildEdges())
            getNodesRecurse(childLineage, networkNodeSet);
    }

    /**
     * @return set of leaf nodes in network
     */
    public Set<NetworkNode> getLeafNodes() {
        return getNodes().stream().filter(NetworkNode::isLeaf).collect(Collectors.toSet());
    }

    /**
     * @return set of internal nodes in network
     */
    public Set<NetworkNode> getInternalNodes() {
        return getNodes().stream().filter(n -> !n.isLeaf()).collect(Collectors.toSet());
    }

    /**
     * @return set of edge objects comprising network
     */
    public Set<NetworkEdge> getEdges() {
        Set<NetworkEdge> networkEdgeSet = new HashSet<>();

        getEdgesRecurse(rootEdge, networkEdgeSet);

        return networkEdgeSet;
    }

    /**
     * @return number of segments represented on network
     */
    public int getSegmentCount() {
        if (segmentCount == null)
            segmentCount = getLeafNodes().iterator().next().getParentEdges().get(0).hasSegments.cardinality();

        return segmentCount;
    }

    private void getEdgesRecurse(NetworkEdge edge, Set<NetworkEdge> networkEdgeSet) {

        if (networkEdgeSet.contains(edge))
            return;

        networkEdgeSet.add(edge);
        for (NetworkEdge childEdge : edge.childNode.getChildEdges())
            getEdgesRecurse(childEdge, networkEdgeSet);
    }

    /**
     * @return Extended Newick representation of network
     */
    public String getExtendedNewick() {
        return getExtendedNewick(rootEdge, new ArrayList<>(), false) + ";";
    }

    /**
     * @return Extended Newick representation of network, with
     *         segment presence annotation.
     */
    public String getExtendedNewickVerbose() {
        return getExtendedNewick(rootEdge, new ArrayList<>(), true) + ";";
    }

    private String getExtendedNewick(NetworkEdge currentEdge, List<NetworkNode> seenReassortmentNodes,
                                     boolean verbose) {
        StringBuilder result = new StringBuilder();

        boolean traverse = true;
        int hybridID = -1;
        if (currentEdge.childNode.isReassortment()) {
            hybridID = seenReassortmentNodes.indexOf(currentEdge.childNode);

            if (hybridID<0) {
                traverse = false;
                seenReassortmentNodes.add(currentEdge.childNode);
                hybridID = seenReassortmentNodes.size()-1;
            }
        }

        if (traverse && !currentEdge.childNode.isLeaf()) {
            result.append("(");

            boolean isFirst = true;
            for (NetworkEdge childEdge : currentEdge.childNode.getChildEdges()) {
                if (isFirst)
                    isFirst = false;
                else
                    result.append(",");

                result.append(getExtendedNewick(childEdge, seenReassortmentNodes, verbose));
            }

            result.append(")");
        }

        if (currentEdge.childNode.getTaxonLabel() != null)
            result.append(currentEdge.childNode.getTaxonLabel());

        if (hybridID>=0) {
            result.append("#H").append(hybridID);
        }

        result.append("[&");
        result.append("segments=").append(currentEdge.hasSegments);
        if (verbose) {
            for (int segIdx=0; segIdx<getSegmentCount(); segIdx++) {
                result.append(",seg").append(segIdx).append("=")
                        .append(currentEdge.hasSegments.get(segIdx));
            }
        }
        result.append(",segsCarried=").append(currentEdge.hasSegments.cardinality());
        result.append("]");

        if (currentEdge.parentNode != null)
            result.append(":").append(currentEdge.parentNode.getHeight() - currentEdge.childNode.getHeight());
        else
            result.append(":0.0");

        return result.toString();
    }

    public void fromExtendedNewick(String newickStr) {

        CharStream inputStream = CharStreams.fromString(newickStr);
        NetworkLexer lexer = new NetworkLexer(inputStream);
        CommonTokenStream tokenStream = new CommonTokenStream(lexer);
        NetworkParser parser = new NetworkParser(tokenStream);
        ParseTree tree = parser.network();

        NetworkBuilderVisitor builder = new NetworkBuilderVisitor();
        rootEdge = builder.visit(tree);
    }

    /** StateNode implementation: **/

    @Override
	public String toString() {
        return getExtendedNewick();
    }

    @Override
    public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
    }

    @Override
    public StateNode copy() {
        return new Network(copyEdge(rootEdge, new HashMap<>()));
    }

    private NetworkEdge copyEdge(NetworkEdge edge, Map<NetworkNode,NetworkNode> seenNodes) {

        NetworkEdge edgeCopy = new NetworkEdge(null, null, (BitSet)edge.hasSegments.clone());
        NetworkNode childNodeCopy;
        boolean traverse = true;
        if (seenNodes.containsKey(edge.childNode)) {
            childNodeCopy = seenNodes.get(edge.childNode);
            traverse = false;
        } else {
            childNodeCopy = new NetworkNode();
            childNodeCopy.setHeight(edge.childNode.getHeight());
            childNodeCopy.setLabel(edge.childNode.getTaxonLabel());
            seenNodes.put(edge.childNode, childNodeCopy);
        }

        childNodeCopy.addParentEdge(edgeCopy);

        if (traverse) {
            for (NetworkEdge childEdge : edge.childNode.getChildEdges()) {
                NetworkEdge childEdgeCopy = copyEdge(childEdge, seenNodes);
                childNodeCopy.addChildEdge(childEdgeCopy);
            }
        }

        return edgeCopy;
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
        fromExtendedNewick(node.getTextContent().replaceAll("&amp;", "&"));
    }

    @Override
    public int scale(double scale) {
        return 0;
    }

    @Override
    protected void store() {
        storedRootEdge = copyEdge(rootEdge, new HashMap<>());
    }

    @Override
    public void restore() {
        NetworkEdge tmp = storedRootEdge;
        storedRootEdge = rootEdge;
        rootEdge = tmp;
        hasStartedEditing = false;
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

    /** Loggable implementation: **/

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
                NetworkNode childNode = childEdge.childNode;
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
                shiftNodeHeights(maxRTLT, childEdge.childNode, seenNodes);
        }

        @Override
        public NetworkEdge visitNetwork(NetworkParser.NetworkContext ctx) {
            seenHybrids = new HashMap<>();
            edgeLengths = new HashMap<>();

            NetworkEdge rootEdge = visit(ctx.node());

            Set<NetworkNode> seenNodes = new HashSet<>();
            NetworkNode rootNode = rootEdge.childNode;
            rootNode.setHeight(0.0);
            double maxRTLT = getMaxRootToLeafTime(rootNode, seenNodes);

            seenNodes.clear();
            shiftNodeHeights(maxRTLT, rootEdge.childNode, seenNodes);

            return rootEdge;
        }

        private String removeQuotes(String str) {

            String[] quoteChars = {"\"", "'"};

            for (String quoteChar : quoteChars) {
                if (str.startsWith(quoteChar) && str.endsWith(quoteChar) && str.length() >= 2)
                    str = str.substring(1, str.length() - 1);
            }

            return str;
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
                node.setLabel(removeQuotes(ctx.post().label().getText()));

            for (NetworkParser.NodeContext childNodeCtx : ctx.node()) {
                NetworkEdge childEdge = visit(childNodeCtx);
                childEdge.parentNode = node;
                node.addChildEdge(childEdge);
            }

            boolean segmentsProcessed = false;
            BitSet hasSegments = new BitSet();
            if (ctx.post().meta() != null
                    && ctx.post().meta().attrib() != null) {

                for (NetworkParser.AttribContext attribCtx : ctx.post().meta().attrib()) {
                    if (!removeQuotes(attribCtx.attribKey.getText()).equals("segments"))
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
     * Modify segmentTree so that it matches the current network state.
     *
     * @param segmentTree
     */
    public void updateSegmentTree(Tree segmentTree, int segmentIdx) {

    }
}