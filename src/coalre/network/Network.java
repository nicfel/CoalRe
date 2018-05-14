package coalre.network;

import beast.core.StateNode;
import beast.util.XMLParser;
import org.w3c.dom.Node;

import java.io.PrintStream;

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

    public void fromExtendedNewick(String newickString) {

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
        fromExtendedNewick(node.getTextContent().replaceAll("&amp;","&"));
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