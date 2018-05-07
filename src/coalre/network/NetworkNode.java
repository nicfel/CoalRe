package coalre.network;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class NetworkNode {


    /**
     * height of this node.
     */
    protected double height = Double.MAX_VALUE;


    /**
     * list of children of this node *
     * Don't use m_left and m_right directly
     * Use getChildCount() and getChild(x) or getChildren() instead
     */
    List<NetworkNode> children = new ArrayList<>();

    /**
     * parent node in the beast.tree, null if root *
     */
    NetworkNode parent = null;
    NetworkNode secondparent = null;

    Integer reassortmentNumber;
    
    /**
     * node number of the corresponding segment tree *
     */
    Integer[] segNodeNr;
    Integer[] goesFirst;
    Integer[] goesSecond;
    
    Boolean[] hasSegments;

    public NetworkNode() {
    }

    public double getHeight() {
        return height;
    }

    public void setHeight(final double height) {
        this.height = height;
    }

    /**
     * @return length of branch between this node and its parent in the beast.tree
     */
    public final double getLength() {
        if (isRoot()) {
            return 0;
        } else {
            return getParent().height - height;
        }
    }

    /**
     * @return parent node, or null if this is root *
     */
    public NetworkNode getParent() {
        return parent;
    }

    /**
     * Calls setParent(parent, true)
     *
     * @param parent the new parent to be set, must be called from within an operator.
     */
    void setParent(final NetworkNode parent) {
        	this.parent = parent;
    }
    
    public NetworkNode getSecondParent(){
    	return secondparent;
    }
    
    public void setSecondParent(final NetworkNode parent) {
        this.secondparent = parent;
    }
    
    public void setReassortmentNumber(Integer reassortmentNumber){
    	this.reassortmentNumber = reassortmentNumber;
    }
    
    public Integer getReassortmentNumber(){
    	return reassortmentNumber;
    }


    public Integer[] getSegNodes(){
    	return segNodeNr;
    }
    
    public void setSegNodes(Integer[] segNodeNr){
    	//TODO pointer issue?
    	this.segNodeNr = segNodeNr;
    }
    
    
    public Boolean[] getHasSegments(){
    	return hasSegments;
    }
    
    public void setHasSegments(Boolean[] hasSegments){
    	this.hasSegments = new Boolean[hasSegments.length];
    	//TODO there is probably a better way of doing this
    	for (int i = 0; i < hasSegments.length; i++){
    		if (hasSegments[i]==null)
    			this.hasSegments[i] = false;
    		else
    			this.hasSegments[i] = hasSegments[i];
    	}
    }
    
    
    public Integer[] getGoesFirst(){
    	return goesFirst;
    }
    
    public Integer[] getGoesSecond(){
    	return goesSecond;
    }
    
    
    public void setGoesFirst(Integer[] goesFirst){
    	//TODO pointer issue?
    	this.goesFirst = goesFirst;
    }
    
    public void getGoesSecond(Integer[] goesSecond){
    	//TODO pointer issue?
    	this.goesSecond = goesSecond;
    }    
    
    /**
     * @return unmodifiable list of children of this node
     */
    public List<NetworkNode> getChildren() {
        return Collections.unmodifiableList(children);
    }

    /**
     * @return true if current node is root node *
     */
    public boolean isRoot() {
        return parent == null;
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


    public void removeChild(final NetworkNode child) {
        children.remove(child);
    }

    public void addChild(final NetworkNode child) {
        child.setParent(this);
        children.add(child);
    }

    /**
     * @return count number of nodes in beast.tree, starting with current node *
     */
    public int getNodeCount() {
        int nodes = 1;
        for (final NetworkNode child : children) {
            nodes += child.getNodeCount();
        }
        return nodes;
    }

    public int getLeafNodeCount() {
        if (isLeaf()) {
            return 1;
        }
        int nodes = 0;
        for (final NetworkNode child : children) {
            nodes += child.getLeafNodeCount();
        }
        return nodes;
    }

    public int getInternalNodeCount() {
        if (isLeaf()) {
            return 0;
        }
        int nodes = 1;
        for (final NetworkNode child : children) {
            nodes += child.getInternalNodeCount();
        }
        return nodes;
    }

    /**
     * @return the number of children of this node.
     */
    public int getChildCount() {
        return children.size();
    }

    /**
     * This method returns the i'th child, numbering starting from 0.
     *
     * This method is unprotected and will throw an ArrayOutOfBoundsException if provided an index larger than getChildCount() - 1, or smaller than 0.
     *
     * @return the i'th child of this node.
     */
    public NetworkNode getChild(final int childIndex) {
        return children.get(childIndex);
    }
}
