package coalre.network;

import beast.core.Input;
import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;
import org.w3c.dom.Node;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Network extends StateNode {
    final public Input<TaxonSet> taxonSetInput = new Input<>("taxonset",
            "set of taxa that correspond to the leafs in the network");

    protected NetworkNode root;


    public Network() {
    }

    public Network(final NetworkNode rootNode) {
        setRoot(rootNode);
    }

    @Override
    public void initAndValidate() { }

	public NetworkNode getRoot() {
        return root;
    }

    public void setRoot(final NetworkNode root) {
        this.root = root;
    }

    public Set<NetworkNode> getNodes() {
        Set<NetworkNode> nodeSet = new HashSet<>();
        addDescendentsToNodeSet(root, nodeSet);

        return nodeSet;
    }

    private void addDescendentsToNodeSet(NetworkNode node, Set<NetworkNode> nodeSet) {
        if (!nodeSet.contains(node)) {
            nodeSet.add(node);

            for (NetworkNode child : node.getChildren())
                addDescendentsToNodeSet(child, nodeSet);
        }
    }

    /** StateNode implementation: **/

    @Override
	public String toString() {
        return root.toString();
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

    @Override
    public void init(PrintStream out) {

    }

    @Override
    public void close(PrintStream out) {

    }
}