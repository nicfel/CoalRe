package coalre.statistics;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Tree;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

public class NetworkStatsLogger extends BEASTObject implements Loggable {


    public Input<Network> networkInput = new Input<>("network",
            "Network for which to log statistics.",
            Input.Validate.REQUIRED);
    
    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "Segment tree associated with network.",
            new ArrayList<>());


    Network network;
    boolean logObservable = false;

    public NetworkStatsLogger() { }

    @Override
    public void initAndValidate() {
        network = networkInput.get();
        if(segmentTreesInput.get().size()>0)
        	logObservable = true;
    }

    @Override
    public void init(PrintStream out) {

        String prefix = network.getID() == null ? "networkStat." : network.getID() + ".";

        if (logObservable)
	        out.print(prefix + "obsHeight\t" +
	                prefix + "obsTotalLength\t" +
	                prefix + "obsReassortmentNodeCount\t");
        else
        	out.print(prefix + "height\t" +
	                prefix + "totalLength\t" +
	                prefix + "reassortmentNodeCount\t");

    }

    @Override
    public void log(long sample, PrintStream out) {
        if (logObservable){
    		double[] rootHeights = new double[segmentTreesInput.get().size()];
			for (int i = 0; i < segmentTreesInput.get().size(); i++)
				rootHeights[i] = segmentTreesInput.get().get(i).getRoot().getHeight();
			
        	out.print(getTotalHeight(network, rootHeights) + "\t" +
	                getTotalEdgeLength(network, rootHeights) + "\t" +
	                getReassortmentCount(network, rootHeights) + "\t");

        	
        }else{   
        	out.print(getTotalHeight(network) + "\t" +
	                getTotalEdgeLength(network) + "\t" +
	                getReassortmentCount(network) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {

    }

    public static int getReassortmentCount(Network network) {
        return (int)network.getNodes().stream().filter(NetworkNode::isReassortment).count();
    }

    public static double getTotalEdgeLength(Network network) {
        return network.getEdges().stream().filter(e -> !e.isRootEdge()).
                map(NetworkEdge::getLength).reduce((l1, l2) -> l1+l2).get();
    }

    public static double getTotalHeight(Network network) {
        return network.getRootEdge().childNode.getHeight();
    }
    
    public static int getReassortmentCount(Network network, double[] rootHeights) {
    	double maxHeight = 0.0;
    	for (int i = 0; i < rootHeights.length; i++)
    		if (rootHeights[i] > maxHeight)
    			maxHeight = rootHeights[i];
    	
    	final double finalMaxHeight = maxHeight;

        return (int)network.getNodes().stream()
        		.filter(NetworkNode::isReassortment)
        		.filter(n -> n.getHeight() < finalMaxHeight)
        		.count();
    }

    public static double getTotalEdgeLength(Network network, double[] rootHeights) {
    	double totalLength = 0.0;
    	List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());
        List<NetworkEdge> nonRootEdges = networkEdges.stream()
        		.filter(e -> !e.isRootEdge())
                .collect(Collectors.toList());
        // check for each edge if it has at least one segment for which the root hasn't been reched yet
        for (int i = 0; i < nonRootEdges.size(); i++){
        	double childHeight = nonRootEdges.get(i).childNode.getHeight();
    		final BitSet hasSegment = nonRootEdges.get(i).hasSegments;
    		for (int j = 0; j < rootHeights.length; j++){
    			if (hasSegment.get(j) && childHeight < rootHeights[j]){
    				totalLength += nonRootEdges.get(i).parentNode.getHeight() - childHeight;
    				break;
    			}
    		}
        	
        }        	
    	
        return totalLength;
    }

    public static double getTotalHeight(Network network, double[] rootHeights) {
    	double maxHeight = 0.0;
    	for (int i = 0; i < rootHeights.length; i++)
    		if (rootHeights[i] > maxHeight)
    			maxHeight = rootHeights[i];
        return maxHeight;
    }

}
