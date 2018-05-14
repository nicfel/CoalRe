package coalre.network;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class NetworkNode {

    /**
     * Taxon corresponding to this node (usually only used for leaves)
     */
    String taxonLabel;

    /**
     * height of this node.
     */
    protected double height = Double.MAX_VALUE;


    List<NetworkEdge> children = new ArrayList<>();
    List<NetworkEdge> parents = new ArrayList<>();

    Integer reassortmentNumber;

    /**
     * Which segments do the lineages above this node correspond to?
     */
    List<Boolean[]> hasSegments;

    public NetworkNode() {
    }

    public double getHeight() {
        return height;
    }

    public NetworkNode setHeight(final double height) {
        this.height = height;
        return this;
    }

    /**
     * @return parent node, or null if this is root *
     */
    public List<NetworkEdge> getParentEdges() {
        return parents;
    }

    public NetworkNode addParentEdge(NetworkEdge newParentEdge) {
        parents.add(newParentEdge);
        return this;
    }

    public NetworkNode removeParentEdge(NetworkEdge parentEdge) {
        parents.remove(parentEdge);
        return this;
    }

    public List<NetworkEdge> getChildEdges() {
        return children;
    }

    public NetworkNode addChildEdge(NetworkEdge newChildEdge) {
        children.add(newChildEdge);
        return this;
    }

    public NetworkNode removeChildEdge(NetworkEdge childEdge) {
        children.remove(childEdge);
        return this;
    }


    public NetworkNode setReassortmentNumber(Integer reassortmentNumber){
    	this.reassortmentNumber = reassortmentNumber;
    	return this;
    }
    
    public Integer getReassortmentNumber(){
    	return reassortmentNumber;
    }

    /**
     * @return true if current node is root node *
     */
    public boolean isRoot() {
        return parents.isEmpty();
    }

    /**
     * @return true if current node is a leaf node *
     */
    public boolean isLeaf() {
        return children.size() == 0;
        //return getLeft() == null && getRight() == null;
    }
    
    public boolean isReassortment() {
        return reassortmentNumber != null;
    }

    public int getChildCount() {
        return children.size();
    }

    public int getParentCount() {
        return parents.size();
    }

    /**
     * @return the label of the taxon (if any) corresponding to this node
     */
    public String getTaxonLabel() {
        return taxonLabel;
    }

    /**
     * Sets the label of the taxon corresponding to this node.
     *
     * @param taxonLabel the new taxon label
     */
    public void setTaxonLabel(String taxonLabel) {
        this.taxonLabel = taxonLabel;
    }
}
