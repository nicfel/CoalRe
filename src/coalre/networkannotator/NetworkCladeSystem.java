package coalre.networkannotator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.evolution.tree.Node;
import beast.math.statistic.DiscreteStatistics;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;


/**
 * extracted from TreeAnnotator
 */
public class NetworkCladeSystem {

    protected List<Map<BitSet, Clade>> cladeMap = new ArrayList<>();
    public List<String> leafNodeMap;


    public NetworkCladeSystem() { 
    }
    
    public void setLeafLabels(List<NetworkNode> leafNodes, int nrSegments){
    	for (int i = 0; i < nrSegments; i++)
    		cladeMap.add(new HashMap<>());

    	leafNodeMap = new ArrayList<String>();
    	for (NetworkNode leaf : leafNodes)
    		leafNodeMap.add(leaf.getTaxonLabel());    	
    }


    /**
     * adds all the clades in the tree
     */
    public void add(Network network, boolean includeTips) {
        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).
    	for (int i = 0; i < network.getSegmentCount(); i++){
    		addClades(network.getRootEdge().childNode, includeTips, i);
    	}
    }    

    private BitSet addClades(NetworkNode node, boolean includeTips, int segment) {
        BitSet bits = new BitSet();

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
        	if (node.getParentEdges().get(0).hasSegments.get(segment))
        		bits.set(index);


        } else {

        	// get the children of that node
            List<NetworkEdge> childEdges = node.getChildEdges();
            // add all children to the bitset
            for (NetworkEdge childEdge : childEdges)
            	if (childEdge.hasSegments.get(segment))
            		bits.or(addClades(childEdge.childNode, includeTips, segment));
            
            if (node.isCoalescence() && 
            		node.getChildEdges().get(0).hasSegments.get(segment) &&
            			node.getChildEdges().get(1).hasSegments.get(segment)){
            	addCoalescentClade(bits, segment);
            }
        }

        return bits;
    }

    private void addCoalescentClade(BitSet bits, int segment) {
        Clade clade = cladeMap.get(segment).get(bits);
        if (clade == null) {
            clade = new Clade(bits);
            cladeMap.get(segment).put(bits, clade);
        }
        clade.setCount(clade.getCount() + 1);
    }
    
    
    public void collectAttributes(Network network, Set<String> attributeNames) {
    	for (int i = 0; i < network.getSegmentCount(); i++){
    		collectAttributes(network.getRootEdge().childNode, attributeNames, i);
    	}
    }

    private BitSet collectAttributes(NetworkNode node, Set<String> attributeNames, int segment) {

        BitSet bits = new BitSet();

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            if (index < 0) {
                throw new IllegalArgumentException("Taxon with height= " + node.getHeight() + ", not found in target tree");
            }
        	if (node.getParentEdges().get(0).hasSegments.get(segment))
        		bits.set(index);

        } else {
           	
        	// get the children of that node
            List<NetworkEdge> childEdges = node.getChildEdges();
            // add all children to the bitset
            for (NetworkEdge childEdge : childEdges)
            	if (childEdge.hasSegments.get(segment))
            		bits.or(collectAttributes(childEdge.childNode, attributeNames, segment));
            
            if (node.isCoalescence() && 
            		node.getChildEdges().get(0).hasSegments.get(segment) &&
            			node.getChildEdges().get(1).hasSegments.get(segment)){
            	collectAttributesForClade(bits, node, attributeNames, segment);
            }
        }

        return bits;
    }

    private void collectAttributesForClade(BitSet bits, NetworkNode node, Set<String> attributeNames, int segment) {
        Clade clade = cladeMap.get(segment).get(bits);
        if (clade != null) {

            if (clade.attributeValues == null) {
                clade.attributeValues = new ArrayList<>();
            }

            int i = 0;
            Object[] values = new Object[attributeNames.size()];
            for (String attributeName : attributeNames) {

                Object value;
                switch (attributeName) {
                    case "height":
                        value = node.getHeight();
                        break;
                    case "length":
                        value = getBranchLength(node);
                        break;
                    default:
                    	throw new IllegalArgumentException("Summary not implemented for values other than height and posterior");
//                        value = node.getMetaData(attributeName);
//                        if (value instanceof String && ((String) value).startsWith("\"")) {
//                            value = ((String) value).replaceAll("\"", "");
//                        }
//                        break;
                }

                values[i] = value;

                i++;
            }
            clade.attributeValues.add(values);

            clade.setCount(clade.getCount() + 1);
        }
    }

    
    
    public void summarizeAttributes(Network network, Set<String> attributeNames, boolean useMean, int nrNetworks) {
    	summarizeAttributes(network.getRootEdge().childNode, attributeNames, useMean, nrNetworks);
    }

    private List<BitSet> summarizeAttributes(NetworkNode node, Set<String> attributeNames, boolean useMean, int nrNetworks) {

        List<BitSet> bits = new ArrayList<BitSet>();
        for (int segment = 0; segment < cladeMap.size(); segment++){
        	bits.add(new BitSet());
        }
        

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            if (index < 0) {
                throw new IllegalArgumentException("Taxon with height= " + node.getHeight() + ", not found in target tree");
            }
            for (int segment = 0; segment < cladeMap.size(); segment++){
            	if (node.getParentEdges().get(0).hasSegments.get(segment))
            		bits.get(segment).set(index);
            }

        } else {
           	
        	// get the children of that node
            List<NetworkEdge> childEdges = node.getChildEdges();
            // add all children to the bitset
            List<BitSet> newbits = new ArrayList<>();
            for (NetworkEdge childEdge : childEdges){
            	newbits = summarizeAttributes(childEdge.childNode, attributeNames, useMean, nrNetworks);

            	for (int segment = 0; segment < cladeMap.size(); segment++){
	            	if (childEdge.hasSegments.get(segment))
	            		bits.get(segment).or(newbits.get(segment));
	            }           	

            } 
            

        }       

        
        if (node.isCoalescence())
        	summarizeAttributesForClade(bits, node, attributeNames, useMean, nrNetworks);
       

        return bits;
    }

    private void summarizeAttributesForClade(List<BitSet> bits, NetworkNode node, 
    		Set<String> attributeNames, boolean useMean, int nrNetworks) {

    	List<NetworkEdge> childEdges = node.getChildEdges();
    	
        for (String attributeName : attributeNames) {
            switch (attributeName) {
            	case "height":
            		node.setMetaData(",segposterior={");
            		
            		int avg_post = 0;
            		int segcount = 0;
            		
            		List<Double> height = new ArrayList<>();

            		for (int segment = 0; segment < cladeMap.size(); segment++){
            			if (segment>0)
				            node.setMetaData(node.getMetaData() + ",");

            			if (childEdges.get(0).hasSegments.get(segment) &&
            					childEdges.get(1).hasSegments.get(segment)){
            				
				    		List<Object[]> rawHeights = cladeMap.get(segment).get(bits.get(segment)).getAttributeValues();
				            double[] heights = new double[rawHeights.size()];
				            for (int i = 0; i < rawHeights.size(); i++)
				            	height.add((double) rawHeights.get(i)[0]);
				            
				            segcount++;
				            avg_post+=rawHeights.size();
				            
				            node.setMetaData(node.getMetaData() + (double)rawHeights.size()/(double)nrNetworks);
//				            if (useMean){
//				                node.setMetaData(node.getMetaData() + ",posterior." + segment + "="
//				                		+DiscreteStatistics.mean(heights));
//				            }else{
//				                node.setHeight(DiscreteStatistics.median(heights));
//				            }
						}else{
				            node.setMetaData(node.getMetaData() + "NA");
						}
					}
            		double avg_pos_val = (double) avg_post/ (double) (nrNetworks*segcount);
            		if (Double.isNaN(avg_pos_val)) avg_pos_val = 0.0;
            		
            		// Convert height to Array
            		double[] heightarray = new double[height.size()];
            		for (int i = 0; i < height.size(); i++)
            			heightarray[i] = height.get(i);
            		
            		if (heightarray.length>0){
						if (useMean){
							node.setHeight(DiscreteStatistics.mean(heightarray));
						}else{
						    node.setHeight(DiscreteStatistics.median(heightarray));
						}
            		}

            		double minHPD,maxHPD;
            		if (heightarray.length>0){
	                    Arrays.sort(heightarray);
	                    minHPD = heightarray[(int)(0.025 * heightarray.length)];
	                    maxHPD = heightarray[(int)(0.975 * heightarray.length)];
            		}else{
            			minHPD = node.getHeight();
            			maxHPD = node.getHeight();
            		}
            		
		            node.setMetaData(",avg_posterior=" + avg_pos_val +
		            		",height_95%_HPD={" + minHPD + "," + maxHPD + "}" + 
		            		node.getMetaData() + "}");

                   
                case "length":
                    break;
                default:
                	throw new IllegalArgumentException("");
            }
        }
    	
    }
    
    private Object getBranchLength(NetworkNode node) {
        if (node.isRoot()) {
            return 0;
        }
        return node.getParent().getHeight() - node.getHeight();
    }

    
    public Map<BitSet, Clade> getCoalescentCladeMap() {
        return coalescentCladeMap;
    }
    
    public Map<BitSet, Clade> getReticulationCladeMap() {
        return reticulationCladeMap;
    }

    public void calculateCladeCredibilities(int nrSegments, int totalTreesUsed) {
    	for (int i = 0; i < nrSegments; i++){
	        for (Clade clade : cladeMap.get(i).values()) {	
	        	
	            if (clade.getCount() > totalTreesUsed) {
	            	for (int j = 0; j < leafNodeMap.size();j++)
	            		if (clade.bits.get(j))
	            			System.err.println(leafNodeMap.get(j));
	
	                throw new AssertionError("clade.getCount=(" + clade.getCount() +
	                        ") should be <= totalTreesUsed = (" + totalTreesUsed + ")");
	            }
	            clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);
	        }
    	}
    }
    

    public double getSumCladeCredibility(NetworkNode node, BitSet bits) {

        double sum = 0.0;

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            bits.set(2*index);
        } else {

            BitSet bits2 = new BitSet();
            for (int i = 0; i < node.getChildCount(); i++) {

            	NetworkNode node1 = node.getChild(i);

                sum += getSumCladeCredibility(node1, bits2);
            }

            for (int i=1; i<bits2.length(); i=i+2) {
                bits2.set(i, false);
            }

            if (node.isFake() && processSA) {
                int index = getTaxonIndex(node.getDirectAncestorChild());
                bits2.set(2 * index + 1);
            }

            sum += getCladeCredibility(bits2);

            if (bits != null) {
                bits.or(bits2);
            }
        }

        return sum;
    }

    public double getLogCladeCredibility(int segment, NetworkNode node, BitSet bits) {

        double logCladeCredibility = 0.0;
        

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            bits.set(index);
        } else {
        	
            BitSet bits2 = new BitSet();
            
        	// get the children of that node
            List<NetworkEdge> childEdges = node.getChildEdges();
            
            // add all children to the bitset
            for (NetworkEdge childEdge : childEdges)
            	if (childEdge.hasSegments.get(segment))
            		logCladeCredibility += getLogCladeCredibility(segment, childEdge.childNode, bits2);
                        
        	if (node.isCoalescence() && 
        			node.getChildEdges().get(0).hasSegments.get(segment) &&
        				node.getChildEdges().get(1).hasSegments.get(segment)){
        			logCladeCredibility += Math.log(getCladeCredibility(segment, bits2));
        	}
            
//            if (logCladeCredibility==Double.NEGATIVE_INFINITY)
            if (bits != null) {
            	bits.or(bits2);
//            	if (node.isCoalescence() && 
//                		node.getChildEdges().get(0).hasSegments.get(segment) &&
//                			node.getChildEdges().get(1).hasSegments.get(segment))
//            		bits.or(bits2);
            }
        }

        return logCladeCredibility;
    }

    private double getCladeCredibility(int segment, BitSet bits) {
        Clade clade = cladeMap.get(segment).get(bits);
        if (clade == null) {
            return 0.0;
        }
        return clade.getCredibility();
    }

    public BitSet removeClades(NetworkNode node, boolean includeTips) {

        BitSet bits = new BitSet();

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            bits.set(2*index);

            if (includeTips) {
                removeClade(bits);
            }

        } else {

            for (int i = 0; i < node.getChildCount(); i++) {

            	NetworkNode node1 = node.getChild(i);

                bits.or(removeClades(node1, includeTips));
            }

            for (int i=1; i<bits.length(); i=i+2) {
                bits.set(i, false);
            }
            if (node.isFake() && processSA) {
                int index = getTaxonIndex(node.getDirectAncestorChild());
                bits.set(2 * index + 1);
            }

            removeClade(bits);
        }

        return bits;
    }

    private void removeClade(BitSet bits) {
        Clade clade = cladeMap.get(bits);
        if (clade != null) {
            clade.setCount(clade.getCount() - 1);
        }

    }

    // Get tree clades as bitSets on target taxa
    // codes is an array of existing BitSet objects, which are reused

//    void getTreeCladeCodes(Network tree, BitSet[] codes) {
//        getTreeCladeCodes(tree.getRootEdge().childNode, codes);
//    }
//
//    int getTreeCladeCodes(NetworkNode node, BitSet[] codes) {
//        final int inode = node.getNr();
//        codes[inode].clear();
//        if (node.isLeaf()) {
//            int index = getTaxonIndex(node);//getTaxonIndex(node);
//            codes[inode].set(index);
//        } else {
//            for (int i = 0; i < node.getChildCount(); i++) {
//                final NetworkNode child = node.getChild(i);
//                final int childIndex = getTreeCladeCodes(child, codes);
//
//                codes[inode].or(codes[childIndex]);
//            }
//        }
//        return inode;
//    }
//    
    /**
     * get the index of a leaf node
     */ 
    private int getTaxonIndex(NetworkNode leaf){
    	return leafNodeMap.indexOf(leaf.getTaxonLabel());
    }

    public class Clade {
        public Clade(BitSet bits) {
            this.bits = bits;
            count = 0;
            credibility = 0.0;
        }

        public int getCount() {
            return count;
        }

        public void setCount(int count) {
            this.count = count;
        }

        public double getCredibility() {
            return credibility;
        }

        public void setCredibility(double credibility) {
            this.credibility = credibility;
        }

        public List<Object[]> getAttributeValues() {
            return attributeValues;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final Clade clade = (Clade) o;

            return !(bits != null ? !bits.equals(clade.bits) : clade.bits != null);

        }

        @Override
        public int hashCode() {
            return (bits != null ? bits.hashCode() : 0);
        }

        @Override
        public String toString() {
            return "clade " + bits.toString() + " #" + count + " count " + getCount();
        }

        int count;
        double credibility;
        BitSet bits;
        List<Object[]> attributeValues = null;
    }

}
