package coalre.network;

import beast.core.StateNode;
import org.w3c.dom.Node;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Network extends StateNode {

    protected NetworkEdge rootEdge;


    public Network() {
    }

    public Network(NetworkEdge rootEdge) {
        setRootEdge(rootEdge);
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

    @Override
    public void init(PrintStream out) {

    }

    @Override
    public void close(PrintStream out) {

    }
}