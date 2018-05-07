package coalre.network;



import java.util.List;

import beast.evolution.alignment.TaxonSet;

public interface NetworkInterface {
    String getID();

    int getLeafNodeCount();
	int getInternalNodeCount();
	int getNodeCount();

	Networknode getRoot();
	Networknode getNode(int i);
	Networknode [] getNodesAsArray();

    List<Networknode> getExternalNodes();
    List<Networknode> getInternalNodes();
    
    TaxonSet getTaxonset();
    
	boolean somethingIsDirty();

    public void getMetaData(Networknode node, Double[] t, String pattern);
    public void setMetaData(Networknode node, Double[] t, String pattern);

    /*
    * Note that leaf nodes are always numbered 0,...,nodeCount-1
    * Internal nodes are numbered higher, but the root has no guaranteed 
    * number.
    */


    /**
     * @param node  top of tree/sub tree (null defaults to whole tree)
     * @param nodes array to fill (null will result in creating a new one)
     * @return tree nodes in post-order, children before parents
     */
    default Networknode[] listNodesPostOrder(Networknode node, Networknode[] nodes) {
        if (node == null) {
            node = getRoot();
        }
        if (nodes == null) {
            // overall node count is cached, faster
            final int n = node == getRoot() ? getNodeCount() : node.getNodeCount();
            nodes = new Networknode[n];
        }
        getNodesPostOrder(node, nodes, 0);
        return nodes;
    }

    static int
    getNodesPostOrder(final Networknode node, final Networknode[] nodes, int pos) {
        //node.m_tree = this;
        for (final Networknode child : node.children) {
            pos = getNodesPostOrder(child, nodes, pos);
        }
        nodes[pos] = node;
        return pos + 1;
    }
}
