package coalre.networkannotator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import beast.base.util.DiscreteStatistics;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.networkannotator.NetworkCladeSystem.BitSetArray;
import coalre.networkannotator.NetworkCladeSystem.DummyClade;
import coalre.networkannotator.NetworkCladeSystem.ReassortmentClade;

public class ReassortmentEventCladeSystem extends NetworkCladeSystem {
	
	
    public void add(Network network, boolean includeTips, Set<String> attributeNames) {
        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).   
    	newCoalescentCladeMap = new ArrayList<>();
		newReassortmentCladeMap = new ArrayList<>();
		for (int i = 0; i < network.getSegmentCount(); i++){
			newCoalescentCladeMap.add(new HashMap<>());
			newReassortmentCladeMap.add(new HashMap<>());
			started = false;
			addReassortmentCladesOnly(network.getRootEdge().childNode, includeTips, i, attributeNames);
		}
		
		// build the reassortment clades with all segments
		buildReassortmentCladeMap(network.getSegmentCount(), attributeNames);	
    }    
	
    private BitSet addReassortmentCladesOnly(NetworkNode node, boolean includeTips, int segment, Set<String> attributeNames) {
    	BitSet bits = new BitSet();
    	
    	if (!started && node.isCoalescence())
    		if (node.getChildEdges().get(0).hasSegments.get(segment) &&
    				node.getChildEdges().get(1).hasSegments.get(segment) &&
    				node.getParentEdges().get(0).hasSegments.get(segment))
    			return null;

    	
    	if (!started && node.isCoalescence())
    		if (node.getChildEdges().get(0).hasSegments.get(segment) &&
    				node.getChildEdges().get(1).hasSegments.get(segment) &&
    				!node.getParentEdges().get(0).hasSegments.get(segment))
    			started = true;
    		
        

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
        	if (node.getParentEdges().get(0).hasSegments.get(segment))
    			bits.set(index);

        } else {

        	// get the children of that node
            List<NetworkEdge> childEdges = node.getChildEdges();
            
            // add all children to the bitset
            for (NetworkEdge childEdge : childEdges){
            	if (started){
		        	if (childEdge.hasSegments.get(segment))
		    			bits.or(addReassortmentCladesOnly(childEdge.childNode, includeTips, segment, attributeNames));
            	}else{
            		addReassortmentCladesOnly(childEdge.childNode, includeTips, segment, attributeNames);
            	}
            	
            }
            
            // if node is coalescent, add the bitset if the coalescent event is observed on a segment tree
            if (node.isCoalescence() && started){
            }else if (started){
            	if (node.getParentEdges().get(0).hasSegments.get(segment))
            		addReassortmentClade(bits, segment, node.getHeight(), false, attributeNames);
            	else if (node.getParentEdges().get(1).hasSegments.get(segment))
            		addReassortmentClade(bits, segment, node.getHeight(), true, attributeNames);
            }              
        }

        return bits;
    }
    
    /**
     * adds new dummy clades to the new clade map
     * @param bits
     * @param segment
     * @param height
     */
    protected void addReassortmentClade(BitSet bits, int segment, Double height, boolean isLeft, Set<String> attributeNames) {
        DummyClade clade = newReassortmentCladeMap.get(segment).get(height);
        if (clade == null) {
            clade = new DummyClade(bits, height, isLeft);
            
            clade.attributeValues = new ArrayList<>();
            
            int i = 0;
            Object[] values = new Object[attributeNames.size()];

            for (String attributeName : attributeNames) {

                Object value;
                switch (attributeName) {
                    case "height":
                        value = height;
                        break;
                    case "length":
//                        value = getBranchLength(node);
//                        break;
                    default:
                    	throw new IllegalArgumentException("Summary not implemented for values other than height and posterior");
                }

                values[i] = value;

                i++;
            }
            clade.attributeValues.add(values);

            
            newReassortmentCladeMap.get(segment).put(height, clade);
        }else{
        	throw new IllegalArgumentException("reassortment clade should never be found");
        }
    }  
    
    protected void buildReassortmentCladeMap(int nrSegments, Set<String> attributeNames) {
    	// get all unique reassortment node heights
    	List<Double> nodeHeights = new ArrayList<>();
    	List<Boolean[]> segmentDirection = new ArrayList<>();
    	for (int i = 0; i < newReassortmentCladeMap.size(); i++){
    		for (Double height : newReassortmentCladeMap.get(i).keySet()){
    			int index = nodeHeights.indexOf(height);
    			if (index==-1){
    				nodeHeights.add(height);
    				Boolean[] newSegs = new Boolean[newReassortmentCladeMap.size()];
    				newSegs[i] = newReassortmentCladeMap.get(i).get(height).isLeft;
    				segmentDirection.add(newSegs);
    			}else{
    				Boolean[] newSegs = segmentDirection.get(index);
    				newSegs[i] = newReassortmentCladeMap.get(i).get(height).isLeft;
    				segmentDirection.set(index, newSegs);
    			}
    		}
    	}    	
    	
    	for (int i = 0; i < nodeHeights.size(); i++){
    		// check if the first segment that is involved in the reassortment event is going left or right
    		Boolean isLeft = false;
    		for (Boolean dir : segmentDirection.get(i)){
    			if (dir !=null && dir==true){
    				isLeft = true;
    				break;
    			}else if (dir !=null &&dir==false){
    				isLeft = false;
    				break;
    			}  
    					
    		}
    		
    		
    		// make a bit set array that is empty if a segment goes left
    		BitSet[] bits = new BitSet[nrSegments*2];
    		for (int j = 0; j < nrSegments; j++){
    			if (segmentDirection.get(i)[j]!=null){
    				if (isLeft && segmentDirection.get(i)[j])
    					bits[2*j] = newReassortmentCladeMap.get(j).get(nodeHeights.get(i)).bits;
    				else if (isLeft && !segmentDirection.get(i)[j])
    					bits[2*j+1] = newReassortmentCladeMap.get(j).get(nodeHeights.get(i)).bits;
    				else if (!isLeft && !segmentDirection.get(i)[j])
    					bits[2*j] = newReassortmentCladeMap.get(j).get(nodeHeights.get(i)).bits;
    				else
    					bits[2*j+1] = newReassortmentCladeMap.get(j).get(nodeHeights.get(i)).bits;

    			}
    		}
    		
    		BitSetArray bitsArray = new BitSetArray(bits);
    		
    		// add the bits to a new reassortmentClade
            ReassortmentClade clade = reassortmentCladeMap.get(bitsArray);
            if (clade == null) {
                clade = new ReassortmentClade(bits);
                
            	if (clade.attributeValues==null)
            		clade.attributeValues = new ArrayList<>();
            		
            	// add the attributes
            	for (int j = 0; j < nrSegments; j++){
            		if (newReassortmentCladeMap.get(j).get(nodeHeights.get(i))!=null){
            			clade.attributeValues.addAll(newReassortmentCladeMap.get(j).get(nodeHeights.get(i)).getAttributeValues());
            		}
            	}                
            	reassortmentCladeMap.put(bitsArray, clade);
                
            }else{
            	if (clade.attributeValues==null)
            		throw new IllegalArgumentException("clade doesn't have attributes values");
            		
            	// add the attributes
            	for (int j = 0; j < nrSegments; j++){
            		if (newReassortmentCladeMap.get(j).get(nodeHeights.get(i))!=null){
            			clade.attributeValues.addAll(newReassortmentCladeMap.get(j).get(nodeHeights.get(i)).getAttributeValues());
                	}
            	}
            	



            }
            clade.setCount(clade.getCount() + 1);
    	}
    	
    }   

        
    @Override
    public void calculateCladeCredibilities(int totalTreesUsed) {
        for (ReassortmentClade clade : reassortmentCladeMap.values()) {	        	
            if (clade.getCount() > totalTreesUsed) {
               System.err.println("clade.getCount=(" + clade.getCount() +
                        ") should be <= totalTreesUsed = (" + totalTreesUsed + ") (can happen for reassortment clades)");
            }
            clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);
        }    	
    }
    
    
    public void printToFile(File outFile) throws FileNotFoundException{
        try (PrintStream ps = new PrintStream(outFile)) {
        	// print the header
        	ps.println("segment.left\tsegment.right\tposterior\theight.mean\theight.median\theight.lower\theight.lower\tdescendents");
        	
        	for (BitSetArray key : reassortmentCladeMap.keySet()){
        		// get which segments went left and which went right
        		BitSet[] bits = reassortmentCladeMap.get(key).bits;
        		boolean[] left = new boolean[bits.length/2];
        		boolean[] right = new boolean[bits.length/2];
        		
        		for (int i = 0; i < bits.length/2; i++){
        			if (bits[2*i]!=null)
        				left[i] = true;
        			if (bits[2*i+1]!=null)	
        				right[i] = true;
        		}
        		ps.print(Arrays.toString(left) + "\t" + Arrays.toString(right) + "\t");
        		ps.print(reassortmentCladeMap.get(key).getCredibility() +"\t");
        		
        		
        		// get the heights
        		List<Double> height = new ArrayList<>();        		
	    		List<Object[]> rawHeights = reassortmentCladeMap.get(key).getAttributeValues();
	            for (int i = 0; i < rawHeights.size(); i++)
	            	height.add((double) rawHeights.get(i)[0]);

        		// Convert height to Array
        		double[] heightarray = new double[height.size()];
        		for (int i = 0; i < height.size(); i++)
        			heightarray[i] = height.get(i);

        		// compute the heights
        		double minHPD,maxHPD;
                Arrays.sort(heightarray);
                minHPD = heightarray[(int)(0.025 * heightarray.length)];
                maxHPD = heightarray[(int)(0.975 * heightarray.length)];

        		
                ps.print(DiscreteStatistics.mean(heightarray) + "\t");
                ps.print(DiscreteStatistics.median(heightarray) + "\t");
                ps.print(minHPD + "\t");
                ps.print(maxHPD + "\n");
        		
        	}
            ps.close();
        }

    }
    
}
