package coalre.statistics;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Tree;
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
	                getObsReassortmentCount(network) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {

    }
    
    public static int getObsReassortmentCount(Network network) {
        return (int)network.getNodes().stream()
        		.filter(e -> e.isReassortment())
        		.filter(e -> e.getParentEdges().get(0).hasSegments.cardinality()>0)
        		.filter(e -> e.getParentEdges().get(1).hasSegments.cardinality()>0)
        		.count();
    }


    public static int getReassortmentCount(Network network) {
        return (int)network.getNodes().stream().filter(NetworkNode::isReassortment).count();
    }

    public static double getTotalEdgeLength(Network network) {
        return network.getEdges().stream().filter(e -> !e.isRootEdge())
        		.filter(e -> e.hasSegments.cardinality() > 0)
                .map(NetworkEdge::getLength).reduce((l1, l2) -> l1+l2).get();
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

    public static double getLociMRCA(Network network){
        double maxHeight = 0.0;
        for (int i = 0; i < network.getSegmentCount(); i++){
            double height = getHeightSegmentsRoot(network.getRootEdge(), i);
            if (height > maxHeight)
                maxHeight = height;
        }
        return maxHeight;
    }
    
    public static double getSecondHighestLociMRCA(Network network){
        List<Double> rootHeights = new ArrayList<>();
        for (int i = 0; i < network.getSegmentCount(); i++){
        	rootHeights.add(getHeightSegmentsRoot(network.getRootEdge(), i));
        }
        // sort
        rootHeights.sort((h1, h2) -> Double.compare(h2, h1));
        return rootHeights.size() > 1 ? rootHeights.get(1) : rootHeights.get(0);
    }


    static double getHeightSegmentsRoot(NetworkEdge edge, int segment){
        NetworkNode n = edge.childNode;
        if (n.isCoalescence()){
            if (n.getChildEdges().get(0).hasSegments.get(segment)
                    && n.getChildEdges().get(1).hasSegments.get(segment)) {
                return n.getHeight();
            }else if (n.getChildEdges().get(0).hasSegments.get(segment)){
                return getHeightSegmentsRoot(n.getChildEdges().get(0), segment);
            }else{
                return getHeightSegmentsRoot(n.getChildEdges().get(1), segment);
            }
        }else{
            return getHeightSegmentsRoot(n.getChildEdges().get(0), segment);
        }
    }

}
