package coalre.network;

import beast.base.inference.StateNode;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import coalre.network.parser.NetworkBaseVisitor;
import coalre.network.parser.NetworkLexer;
import coalre.network.parser.NetworkParser;
import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTree;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

public class Network extends StateNode {

    protected NetworkEdge rootEdge;

    protected NetworkEdge storedRootEdge;

    public Integer segmentCount = null;
    
    public String[] segmentNames;
    public String baseName;
 
    
    
    public Network() {
    }

    public Network(NetworkEdge rootEdge) {
        this.rootEdge = rootEdge;
    }

    public Network(String newickString) {
        fromExtendedNewick(newickString);
    }

    public Network(String newickString, TaxonSet taxonSet) {
        fromExtendedNewick(newickString);

        for (NetworkNode leafNode : getLeafNodes())
            leafNode.setTaxonIndex(taxonSet.getTaxonIndex(leafNode.getTaxonLabel()));
    }
    

    @Override
    public void initAndValidate() {}

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

//        getEdgesRecurseW(rootEdge, networkEdgeSet);

//        int networkEdgeCount = networkEdgeSet.size();
//        
//        networkEdgeSet = new HashSet<>();
        // Likely faster implementation as it avoids all the contains calls
        getEdgesRecurseNew(rootEdge, networkEdgeSet);
        
//        if (networkEdgeSet.size() != networkEdgeCount)
//            throw new RuntimeException("Network edge count mismatch: " + networkEdgeSet.size() + " != " + networkEdgeCount);
        
        return networkEdgeSet;
    }
    private void getEdgesRecurseW(NetworkEdge edge, Set<NetworkEdge> networkEdgeSet) {
    	getEdgesRecurse(edge, networkEdgeSet);
    }

    private void getEdgesRecurseNewW(NetworkEdge edge, Set<NetworkEdge> networkEdgeSet) {
    	getEdgesRecurseNew(edge, networkEdgeSet);
    }

    
    /**
     * @return number of segments represented on network
     */
    public int getSegmentCount() {
        if (segmentCount == null) {
        	// allow for different number of segments on different leaves
        	BitSet totsegment = new BitSet();
        	for (NetworkNode n : getLeafNodes()) {
        		totsegment.or(n.getParentEdges().get(0).hasSegments);
        	}   
            segmentCount = totsegment.cardinality();
        }

        return segmentCount;
    }

    private void getEdgesRecurse(NetworkEdge edge, Set<NetworkEdge> networkEdgeSet) {

        if (networkEdgeSet.contains(edge))
            return;

        networkEdgeSet.add(edge);
        for (NetworkEdge childEdge : edge.childNode.getChildEdges())
            getEdgesRecurse(childEdge, networkEdgeSet);
    }

    
    private void getEdgesRecurseNew(NetworkEdge edge, Set<NetworkEdge> networkEdgeSet) {

	  	networkEdgeSet.add(edge);
	  	
	  	if (edge.childNode.isReassortment() && edge.equals(edge.childNode.getParentEdges().get(1)))
			return;    		
  		
	    for (NetworkEdge childEdge : edge.childNode.getChildEdges())
	    	getEdgesRecurseNew(childEdge, networkEdgeSet);
    }

    
    
    /**
     * @return Extended Newick representation of network
     */
    public String getExtendedNewick() {
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), false, Integer.MAX_VALUE) + ";";
    }
    
    public String getExtendedNewick(int followSegment) {
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), false, followSegment) + ";";
    }


    /**
     * @return Extended Newick representation of network, with
     *         segment presence annotation.
     */
    public String getExtendedNewickVerbose() {
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), true, Integer.MAX_VALUE) + ";";
    }
    
    public String getExtendedNewickVerbose(int followSegment) {
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), true, followSegment) + ";";
    }


    private String getExtendedNewick(NetworkEdge currentEdge, List<NetworkNode> seenReassortmentNodes, 
    		List<Boolean> isTraverseEdge, boolean verbose, int followSegment) {
        StringBuilder result = new StringBuilder();

        boolean traverse = true;
        int hybridID = -1;
        boolean printMetaData = true;
        if (currentEdge.childNode.isReassortment()) {
        	
            hybridID = seenReassortmentNodes.indexOf(currentEdge.childNode);
            if (hybridID<0) {
	        	List<NetworkEdge> parentEdges = currentEdge.childNode.getParentEdges();
	        	
	        	// get the other edge
	        	NetworkEdge otherEdge;
	        	if (parentEdges.get(0).equals(currentEdge))
	        		otherEdge = parentEdges.get(1);
	        	else
	        		otherEdge = parentEdges.get(0);
	        	
	        	// check which edge is the main edge
	        	if (otherEdge.hasSegments.get(followSegment)){
	                traverse = false;
	                seenReassortmentNodes.add(currentEdge.childNode);
	                isTraverseEdge.add(true);
	                hybridID = seenReassortmentNodes.size()-1;
	        	}else if(currentEdge.hasSegments.get(followSegment)){
	                seenReassortmentNodes.add(otherEdge.childNode);
	                isTraverseEdge.add(false);
	                hybridID = seenReassortmentNodes.size()-1;	        		
	        	}else if (currentEdge.hasSegments.cardinality()<otherEdge.hasSegments.cardinality()){
	                traverse = false;
	                seenReassortmentNodes.add(currentEdge.childNode);
	                isTraverseEdge.add(true);
	                hybridID = seenReassortmentNodes.size()-1;
	        	}else{
	                seenReassortmentNodes.add(otherEdge.childNode);
	                isTraverseEdge.add(false);
	                hybridID = seenReassortmentNodes.size()-1;	        		
	        	}
            }else{
            	traverse = isTraverseEdge.get(hybridID);
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

                result.append(getExtendedNewick(childEdge, seenReassortmentNodes, isTraverseEdge, verbose, followSegment));
            }

            result.append(")");
        }

        if (currentEdge.childNode.getTaxonLabel() != null)
            result.append(currentEdge.childNode.getTaxonLabel());

        if (hybridID>=0) {
            result.append("#H").append(hybridID);
        }

        result.append("[&");
//        if (segmentNames!=null) {
//    		result.append("segments={");
//    		String str = "rem";
//            for (int segIdx=0; segIdx<getSegmentCount(); segIdx++) {
//            	if (currentEdge.hasSegments.get(segIdx)) {
//            		str = str+ ",\"" + segmentNames[segIdx].replace(",", "") + "\"";
//            	}
//            }
//        	result.append(str.replace("rem,", ""));
//        	result.append("}");
//        }else {
        	result.append("segments=").append(currentEdge.hasSegments);
//        }
        if (verbose) {
            for (int segIdx=0; segIdx<getSegmentCount(); segIdx++) {
//                if (segmentNames!=null) {
//	                result.append(",").append(segmentNames[segIdx]).append("=")
//	                        .append(currentEdge.hasSegments.get(segIdx));
//                }else {
	                result.append(",seg").append(segIdx).append("=")
                    	.append(currentEdge.hasSegments.get(segIdx));

//                }
            }
        }
        result.append(",segsCarried=").append(currentEdge.hasSegments.cardinality());
        if (currentEdge.childNode.getTypeLabel() != null) 
        		result.append(",state=").append(currentEdge.childNode.getTypeLabel());
        
        if (currentEdge.childNode.metaDataString != null) 
			result.append(currentEdge.childNode.getMetaData());

        result.append("]");

        if (currentEdge.parentNode != null)
            result.append(":").append(currentEdge.parentNode.getHeight() - currentEdge.childNode.getHeight());
        else
            result.append(":0.0");

        return result.toString();
    }

    public void fromExtendedNewick(String newickStr) {
    	
    	if (newickStr.split(";").length==2) {
    		String segmentString = newickStr.split(";")[0];
    		newickStr = newickStr.split(";")[1] + ";";
    		
    		String[] tmp = segmentString.replace("[", "").replace("]", "").replace(" ", "").split(",");
    		segmentNames = new String[tmp.length-1];
    		baseName = tmp[0];
			for (int i = 1; i < tmp.length;i++) {
				segmentNames[i-1] = tmp[i];
			}
    	}

        CharStream inputStream = CharStreams.fromString(newickStr);
        NetworkLexer lexer = new NetworkLexer(inputStream);
        CommonTokenStream tokenStream = new CommonTokenStream(lexer);
        NetworkParser parser = new NetworkParser(tokenStream);
        ParseTree tree = parser.network();

        NetworkBuilderVisitor builder = new NetworkBuilderVisitor();
        rootEdge = builder.visit(tree);

        List<NetworkNode> leafNodes = new ArrayList<>(getLeafNodes());
    }

    /** StateNode implementation: **/

    @Override
	public String toString() {
    	if (segmentNames!=null) {
    		String[] tmp = new String[segmentNames.length+1];
    		tmp[0] = baseName;
    		for (int i = 0; i < segmentNames.length;i++)
    			tmp[i+1] = segmentNames[i];
    		return getExtendedNewick();
    	}else {
    		return getExtendedNewick();
    	}
    }

    @Override
    public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
    }

    @Override
    public StateNode copy() {
        return new Network(rootEdge.getCopy());
    }

    @Override
    public void assignTo(StateNode other) {
        Network otherNetwork = (Network) other;

        otherNetwork.rootEdge = rootEdge;
        otherNetwork.storedRootEdge = null;
        otherNetwork.segmentCount = null;
    }

    @Override
    public void assignFrom(StateNode other) {
        assignFromFragile(other);
        setID(other.getID());
    }

    @Override
    public void assignFromFragile(StateNode other) {
        Network otherNetwork = (Network) other;

        // Save taxon indices.
        Map<String, Integer> taxonToIndexMap = new HashMap<>();
        getLeafNodes().forEach(n -> taxonToIndexMap.put(n.getTaxonLabel(), n.getTaxonIndex()));

        rootEdge = otherNetwork.rootEdge;
        storedRootEdge = null;
        segmentCount = null;

        // Restore taxon indices
        getLeafNodes().forEach(n -> n.setTaxonIndex(taxonToIndexMap.get(n.getTaxonLabel())));
    }

    @Override
    public void fromXML(org.w3c.dom.Node node) {
        fromExtendedNewick(node.getTextContent().replaceAll("&amp;", "&"));
    }

    @Override
    public int scale(double scale) {
        return 0;
    }

    @Override
    protected void store() {
        storedRootEdge = rootEdge.getCopy();
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
        out.println("Begin trees;");
    }

    @Override
    public void close(PrintStream out) {
        out.println("End;");
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
                node.setTaxonLabel(removeQuotes(ctx.post().label().getText()));

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
                    
//                    if (segmentNames!=null) {
//	                    for (NetworkParser.AttribValueContext attribValueCtx : attribCtx.attribValue().vector().attribValue()) {
//	                    	for (int i = 0; i < segmentNames.length;i++) {
//	                    		if (segmentNames[i].contentEquals(attribValueCtx.getText()))
//	                    			hasSegments.set(i);
//	                    	}
//	                    }
//                    }else {
	                    for (NetworkParser.AttribValueContext attribValueCtx : attribCtx.attribValue().vector().attribValue()) {
	                    	hasSegments.set(Integer.valueOf(attribValueCtx.getText()));
	                    }
//                    }

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

    /** Segment tree updating: **/

    /**
     * Modify segmentTree so that it matches the current network state.
     *
     * @param segmentTree
     */
    public void updateSegmentTree(Tree segmentTree, int segmentIdx) {

        // Identify segment tree clades present in network

        Set<BitSet> cladesInNetwork = new HashSet<>();
        Map<BitSet, NetworkNode> networkCladeNodes = new HashMap<>();
        getSegTreeCladesFromNetwork(getRootEdge().childNode, segmentIdx,
                cladesInNetwork, networkCladeNodes);

        // Identify segment tree clades in existing Tree

        Set<BitSet> cladesInTree = new HashSet<>();
        Map<BitSet, Node> cladeNodes = new HashMap<>();
        getSegTreeClades(segmentTree.getRoot(), cladesInTree, cladeNodes);

        // Identify which clades need to be removed in order to update the
        // Tree to match the segment tree embedded in the network.

        Set<BitSet> removedClades = new HashSet<>(cladesInTree);
        removedClades.removeAll(cladesInNetwork);

        // Remove clades from Tree, placing node objects in a bin

        List<Node> nodeBin = new ArrayList<>();
        for (BitSet clade : removedClades) {
            Node node = cladeNodes.get(clade);
            nodeBin.add(node);

            if (!node.isRoot()) {
                node.getParent().removeChild(node);
                node.setParent(null);
            }

            List<Node> children = new ArrayList<>(node.getChildren());
            for (Node child : children) {
                node.removeChild(child);
                child.setParent(null);
            }
        }

        for (BitSet clade : removedClades)
            cladeNodes.remove(clade);


        // Traverse segment tree in network, adding new nodes from
        // bin as required

        BitSet rootClade = buildSegmentTree(getRootEdge().childNode, segmentIdx, cladeNodes, nodeBin);

        // Keep root up to date

        segmentTree.setRoot(cladeNodes.get(rootClade));
    }

    /**
     * Identify which clades are present in the network below edge for the given segment.
     *
     * @param networkNode       network node below which clades are identified
     * @param segmentIdx        index of segment to consider
     * @param clades            empty set to populate with discovered clades
     * @param cladeNetworkNodes empty map to populate with links from caldes to network nodes.
     * @return clade corresponding to network node
     */
    BitSet getSegTreeCladesFromNetwork(NetworkNode networkNode, int segmentIdx,
                                       Set<BitSet> clades,
                                       Map<BitSet, NetworkNode> cladeNetworkNodes) {

        BitSet thisClade = new BitSet();

        int childrenWithSeg = 0;
        for (NetworkEdge childEdge : networkNode.getChildEdges()) {
            if (childEdge.hasSegments.get(segmentIdx)) {
                thisClade.or(getSegTreeCladesFromNetwork(childEdge.childNode, segmentIdx, clades, cladeNetworkNodes));
                childrenWithSeg += 1;
            }
        }

        if (childrenWithSeg == 0) {
            // Leaf node

            thisClade.set(networkNode.getTaxonIndex());
        }

        if (childrenWithSeg != 1) {
            // Node in segment tree

            clades.add(thisClade);
            cladeNetworkNodes.put(thisClade, networkNode);
        }

        return thisClade;
    }

    /**
     * Cached arrays containing mapping between tree node numbers and
     * network node numbers.
     */
    Map<Tree, int[]> nodeNrArrays = new HashMap<>();

    /**
     * Retrieve an array representing a mapping from segment tree node
     * numbers and network node numbers.  The array has the property
     * that network_node_number = array[tree_node_number].
     *
     * @param segmentTree segment tree
     * @return node number array
     */
    int[] getTreeNodeNumberArray(Tree segmentTree) {
        if (nodeNrArrays.keySet().contains(segmentTree))
            return nodeNrArrays.get(segmentTree);

        int[] nodeNumberMap = new int[segmentTree.getLeafNodeCount()];

        for (int treeNodeNr=0; treeNodeNr<nodeNumberMap.length; treeNodeNr++) {
            Node treeNode = segmentTree.getNode(treeNodeNr);
            String taxonID = segmentTree.getTaxonId(treeNode);

            Optional<NetworkNode> maybeNetworkNode =
                    getLeafNodes().stream()
                            .filter(l -> l.getTaxonLabel().equals(taxonID))
                            .findFirst();

            if (maybeNetworkNode.isPresent())
                nodeNumberMap[treeNodeNr] = maybeNetworkNode.get().getTaxonIndex();
            else
                throw new IllegalArgumentException("Segment tree "
                        + segmentTree.getID() + " contains taxon not present in network.");
        }

        nodeNrArrays.put(segmentTree, nodeNumberMap);
        return nodeNumberMap;
    }

    /**
     * Retrieve network node number for a given segment tree node
     * @param treeNode segment tree node
     * @return corresponding network node number
     */
    int getNetworkNodeNr(Node treeNode) {
        return treeNode.getTree().getTaxonset() == null
                ? treeNode.getNr()
                : getTreeNodeNumberArray(treeNode.getTree())[treeNode.getNr()];
    }

    /**
     * Identify which clades are present in the segment tree below node
     *
     * @param node          node below which clades are identified
     * @param clades        empty set to populate with discovered clades
     * @param cladeNodes    empty map to populate with links from clades to tree nodes
     * @return clade corresponding to segment tree node
     */
    BitSet getSegTreeClades(Node node,
                            Set<BitSet> clades,
                            Map<BitSet, Node> cladeNodes) {

        BitSet thisClade = new BitSet();

        if (node.isLeaf()) {
            thisClade.set(getNetworkNodeNr(node));

        } else {
            for (Node child : node.getChildren()) {
                thisClade.or(getSegTreeClades(child, clades, cladeNodes));
            }
        }

        clades.add(thisClade);
        cladeNodes.put(thisClade, node);

        return thisClade;
    }

    /**
     * Update segment tree below networkNode to match tree implied by network.
     *
     * @param networkNode node below which segment tree is updated
     * @param segmentIdx  index of segment tree
     * @param cladeNodes  map containing links from existing tree clades to tree nodes
     * @param nodeBin     list containing node objects available for new tree clades
     * @return clade corresponding to networkNode
     */
    BitSet buildSegmentTree(NetworkNode networkNode, int segmentIdx,
                            Map<BitSet,Node> cladeNodes, List<Node> nodeBin) {

        BitSet thisClade = new BitSet();
        List<BitSet> childClades = new ArrayList<>();

        for (NetworkEdge childEdge : networkNode.getChildEdges()) {
            if (childEdge.hasSegments.get(segmentIdx)) {
                BitSet childClade = buildSegmentTree(childEdge.childNode, segmentIdx, cladeNodes, nodeBin);
                thisClade.or(childClade);
                childClades.add(childClade);
            }
        }
        
        if (childClades.size() == 1)
        	return thisClade;

        
        if (childClades.isEmpty()) {
            // Leaf node

            thisClade.set(networkNode.getTaxonIndex());

        } else {
        	
            // Internal node

            // Retrieve/create tree node object

            Node treeNode;
            if (!cladeNodes.containsKey(thisClade)) {
                treeNode = nodeBin.get(nodeBin.size() - 1);
                nodeBin.remove(nodeBin.size() - 1);
                cladeNodes.put(thisClade, treeNode);
            } else {
                treeNode = cladeNodes.get(thisClade);
            }

            // Update children (is there a more efficient way to do this??)

            Set<Node> trueChildNodes = new HashSet<>();
            for (BitSet childClade : childClades) {
                Node childNode = cladeNodes.get(childClade);
                if (!treeNode.getChildren().contains(childNode))
                    treeNode.addChild(childNode);

                trueChildNodes.add(childNode);
            }

            List<Node> currentChildNodes = new ArrayList<>(treeNode.getChildren());
            for (Node childNode : currentChildNodes) {
                if (!trueChildNodes.contains(childNode))
                    treeNode.removeChild(childNode);
            }

        }        
        
        // keep track of the node number on the networkNode
        if (networkNode.segmentIndices==null) {
        	networkNode.segmentIndices = new int[getSegmentCount()];
        	// initialise to -1
        	Arrays.fill(networkNode.segmentIndices, -1);
        }
        networkNode.setSegmentIndex(segmentIdx, cladeNodes.get(thisClade).getNr());
        
        if (cladeNodes.get(thisClade).getHeight() != networkNode.getHeight())
            cladeNodes.get(thisClade).setHeight(networkNode.getHeight());
        
        return thisClade;
    }    

    /**
     * Will be used in the next version of BEAST to prevent date trait cloning
     * from breaking the BEAuti model.
     */
    public boolean notCloneable() {
        return true;
    }
}