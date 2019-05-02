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
    protected List<Map<Double, DummyClade>> newReassortmentCladeMap = new ArrayList<>();
    protected Map<BitSetArray, ReassortmentClade> reassortmentCladeMap = new HashMap<>();
    public List<String> leafNodeMap;

    public boolean[] followSegmentAlready;
    private boolean started;
    public NetworkCladeSystem() { 
    }
    
    /**
     * adds all leaf labels to a list
     * @param leafNodes
     * @param nrSegments
     */
    public void setLeafLabels(List<NetworkNode> leafNodes, int nrSegments){
    	for (int i = 0; i < nrSegments; i++){
    		cladeMap.add(new HashMap<>());    		
    	}
    	reassortmentCladeMap = new HashMap<>();

    	
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
		newReassortmentCladeMap = new ArrayList<>();
		for (int i = 0; i < network.getSegmentCount(); i++){
			newReassortmentCladeMap.add(new HashMap<>());
			started = false;
			addClades(network.getRootEdge().childNode, includeTips, i);
		}
		
		// build the reassortment clades with all segments
		buildReassortmentCladeMap(network.getSegmentCount());		
    }    

    private BitSet addClades(NetworkNode node, boolean includeTips, int segment) {
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
		    			bits.or(addClades(childEdge.childNode, includeTips, segment));
            	}else{
            		addClades(childEdge.childNode, includeTips, segment);
            	}
            	
            }
            
            // if node is coalescent, add the bitset if the coalescent event is observed on a segment tree
            if (node.isCoalescence() && started){
        		if (node.getChildEdges().get(0).hasSegments.get(segment) &&
        			node.getChildEdges().get(1).hasSegments.get(segment)){
        			addCoalescentClade(bits, segment);
        		}    			
            }else if (started){
            	if (node.getParentEdges().get(0).hasSegments.get(segment))
            		addReassortmentClade(bits, segment, node.getHeight(), false);
            	else if (node.getParentEdges().get(1).hasSegments.get(segment))
            		addReassortmentClade(bits, segment, node.getHeight(), true);
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
    
    /**
     * adds new dummy clades to the new clade map
     * @param bits
     * @param segment
     * @param height
     */
    private void addReassortmentClade(BitSet bits, int segment, Double height, boolean isLeft) {
        DummyClade clade = newReassortmentCladeMap.get(segment).get(height);
        if (clade == null) {
            clade = new DummyClade(bits, height, isLeft);
            newReassortmentCladeMap.get(segment).put(height, clade);
        }else{
        	throw new IllegalArgumentException("reassortment clade should never be found");
        }
    }   
    
    private void buildReassortmentCladeMap(int nrSegments) {
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
            	reassortmentCladeMap.put(bitsArray, clade);
                
            }              
            	
            clade.setCount(clade.getCount() + 1);


    	}
    	
    }   

    public void collectAttributes(Network network, Set<String> attributeNames) {
		newReassortmentCladeMap = new ArrayList<>();
		for (int i = 0; i < network.getSegmentCount(); i++){
			newReassortmentCladeMap.add(new HashMap<>());
			started = false;
			collectAttributes(network.getRootEdge().childNode, attributeNames, i);		
		}		

		collectAtributesReassortmentCladeMap(network.getSegmentCount(), attributeNames);
    }
    

    private BitSet collectAttributes(NetworkNode node, Set<String> attributeNames, int segment) {
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
            if (index < 0) {
                throw new IllegalArgumentException("Taxon with height= " + node.getHeight() + ", not found in target tree");
            }
        	if (node.getParentEdges().get(0).hasSegments.get(segment))
        			bits.set(index);


        } else {        	
           	
        	// get the children of that node
            List<NetworkEdge> childEdges = node.getChildEdges();
            // add all children to the bitset
            for (NetworkEdge childEdge : childEdges){
            	if (started){
		        	if (childEdge.hasSegments.get(segment))
		    			bits.or(collectAttributes(childEdge.childNode, attributeNames, segment));
            	}else{
            		collectAttributes(childEdge.childNode, attributeNames, segment);
            	}
            }  
            
            if (node.isCoalescence() && started){
        		if(node.getChildEdges().get(0).hasSegments.get(segment) &&
        			node.getChildEdges().get(1).hasSegments.get(segment)){
        			collectAttributesForClade(bits, node, attributeNames, segment);
        		}
            } else if (started) {
            	if (node.getParentEdges().get(0).hasSegments.get(segment))
            		collectAttributesForReassortmentClade(bits, node, attributeNames, segment, node.getHeight(), false);
            	else if (node.getParentEdges().get(1).hasSegments.get(segment))
            		collectAttributesForReassortmentClade(bits, node, attributeNames, segment, node.getHeight(), true);
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
//                    case "length":
//                        value = getBranchLength(node);
//                        break;
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
    

    private void collectAttributesForReassortmentClade(BitSet bits, NetworkNode node,
    		Set<String> attributeNames, int segment, double height,  boolean isLeft) {
    	

        DummyClade clade = newReassortmentCladeMap.get(segment).get(height);
        if (clade == null) {
        	
            clade = new DummyClade(bits, height, isLeft);
        	
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
//                    case "length":
//                        value = getBranchLength(node);
//                        break;
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
            
            newReassortmentCladeMap.get(segment).put(height, clade);
        }else{
        	throw new IllegalArgumentException("reassortment clade should never be found");
        }
    }

    private void collectAtributesReassortmentCladeMap(int nrSegments, Set<String> attributeNames) {
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
            if (clade != null) {
            	if (clade.attributeValues==null)
            		clade.attributeValues = new ArrayList<>();
            		
            	// add the attributes
            	for (int j = 0; j < nrSegments; j++){
            		if (segmentDirection.get(i)[j]!=null){
            			clade.attributeValues.addAll(newReassortmentCladeMap.get(j).get(nodeHeights.get(i)).getAttributeValues());
            		}
            	}
            	clade.setCount(clade.getCount()+1);
            }
        }    	
    }   
    
    public void summarizeAttributes(Network network, Set<String> attributeNames, boolean useMean, int nrNetworks) {
		boolean[] followSegment = new boolean[network.getSegmentCount()];
		for (int i=0;i<network.getSegmentCount();i++)followSegment[i] = false;
		
		followSegmentAlready = Arrays.copyOf(followSegment, followSegment.length);

    	// summarizes all coalescent events
    	summarizeAttributes(network.getRootEdge().childNode, attributeNames, useMean, nrNetworks, network.getSegmentCount(), followSegment);    	
    }
    

    private BitSet[] summarizeAttributes(NetworkNode node, Set<String> attributeNames, boolean useMean, int nrNetworks, int nrSegments, boolean[] followSegment_in) {

        BitSet[] bits = new BitSet[cladeMap.size()];
        for (int i = 0; i < cladeMap.size(); i++) bits[i] = new BitSet();
        
        boolean[] followSegment = Arrays.copyOf(followSegment_in, followSegment_in.length);


        if (node.isLeaf()) {
            int index = getTaxonIndex(node);
            if (index < 0) {
                throw new IllegalArgumentException("Taxon with height= " + node.getHeight() + ", not found in target tree");
            }
            for (int i = 0; i < cladeMap.size(); i++){
            	if (node.getParentEdges().get(0).hasSegments.get(i))
            		bits[i].set(index);
            }
        } else {
        	
            // check if the node is the root of a segment tree, if so, follow the segment
            if (node.isCoalescence()){
	           	 for (int i = 0; i < nrSegments;i++){
	           		if (node.getChildEdges().get(0).hasSegments.get(i) &&
	           			node.getChildEdges().get(1).hasSegments.get(i) &&
	           			!node.getParentEdges().get(0).hasSegments.get(i) &&
	           			!followSegmentAlready[i]){
	           			followSegment[i] = true;
	           			followSegmentAlready[i] = true;
	           		}
	   			}
           }            

           	
        	// get the children of that node
            List<NetworkEdge> childEdges = node.getChildEdges();
            // add all children to the bitset
            for (NetworkEdge childEdge : childEdges){
				boolean[] followSegmentout = Arrays.copyOf(followSegment, followSegment.length);
				for (int i = 0; i < nrSegments;i++){
					if (!childEdge.hasSegments.get(i)){
						followSegmentout[i] = false;
					}
				}	
				BitSet[] newBits = summarizeAttributes(childEdge.childNode, attributeNames, useMean, nrNetworks, nrSegments, followSegmentout);
            	for (int i = 0; i < nrSegments;i++){
	            	if (childEdge.hasSegments.get(i)){
            			bits[i].or(newBits[i]);	   
	            	}
            	}
            }

            if (!node.isReassortment())
            	summarizeAttributesForClade(bits, node, attributeNames, useMean, nrNetworks);   
            else{
            	summarizeAttributesForReassortmentClade(bits, node, attributeNames, useMean, nrNetworks, nrSegments);
            }

        }        
        
        return bits;
    }
    

    private void summarizeAttributesForClade(BitSet[] bits, NetworkNode node, 
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
            				
				    		List<Object[]> rawHeights = cladeMap.get(segment).get(bits[segment]).getAttributeValues();
				            for (int i = 0; i < rawHeights.size(); i++)
				            	height.add((double) rawHeights.get(i)[0]);
				            
				            segcount++;
				            avg_post+=rawHeights.size()-1;
				            
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
            		
		            node.setMetaData(",posterior=" + avg_pos_val +
		            		",height_95%_HPD={" + minHPD + "," + maxHPD + "}" + 
		            		node.getMetaData() + "}");

                   
                case "length":
                    break;
                default:
                	throw new IllegalArgumentException("");
            }
        }
    	
    }
    
    
    private void summarizeAttributesForReassortmentClade(BitSet[] bits, NetworkNode node, 
    		Set<String> attributeNames, boolean useMean, int nrNetworks, int nrSegments) {
    	
    	// check if the first segment goes left or right
    	boolean isLeft = false;
    	for (int i = 0; i < nrSegments; i++){
    		if (node.getParentEdges().get(0).hasSegments.get(i)){
    			isLeft = true;
    			break;
    		}else if (node.getParentEdges().get(1).hasSegments.get(i)){
    			isLeft = false;
    			break;
    		}
    	}
    	
    	BitSet[] bitsarray = new BitSet[nrSegments*2];
    	
    	// Build the bitset array to find the reassortment clade  
    	for (int i = 0; i < nrSegments; i++)
			if (isLeft && node.getParentEdges().get(0).hasSegments.get(i))
				bitsarray[2*i] = bits[i];
			else if (isLeft && node.getParentEdges().get(1).hasSegments.get(i))
				bitsarray[2*i+1] = bits[i];
			else if (!isLeft && node.getParentEdges().get(1).hasSegments.get(i))
				bitsarray[2*i] = bits[i];
			else if (!isLeft && node.getParentEdges().get(0).hasSegments.get(i))
				bitsarray[2*i+1] = bits[i];			

    	BitSetArray keyArray = new BitSetArray(bitsarray);
    	
//    	System.out.println("");
//    	System.out.println(node.getHeight());
//    	System.out.println(Arrays.toString(bitsarray));
//    	System.out.println(nrSegments);
//    	System.out.println(reassortmentCladeMap.get(keyArray).getAttributeValues());
//    	System.exit(0);
    	
    	
        for (String attributeName : attributeNames) {
            switch (attributeName) {
            	case "height":
            		node.setMetaData("");
            		
            		List<Double> height = new ArrayList<>();
            		
		    		List<Object[]> rawHeights = reassortmentCladeMap.get(keyArray).getAttributeValues();
		            for (int i = 0; i < rawHeights.size(); i++)
		            	height.add((double) rawHeights.get(i)[0]);

//		            node.setMetaData(node.getMetaData() + (double)rawHeights.size()/(double)nrNetworks);
            		
            		double posterior = (double) (reassortmentCladeMap.get(keyArray).getCount()-1)/ (double) (nrNetworks);
            		
            		
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
            		
//            		System.out.println(node.getMetaData());
		            node.setMetaData(",posterior=" + posterior +
		            		",height_95%_HPD={" + minHPD + "," + maxHPD + "}" + 
		            		node.getMetaData() + "");

                case "length":
                    break;
                default:
                	throw new IllegalArgumentException("");
            }
        }
    	
    }
    

//    private Object getBranchLength(NetworkNode node) {
//        if (node.isRoot()) {
//            return 0;
//        }
//        return node.getParent().getHeight() - node.getHeight();
//    }
//
//    
//    public Map<BitSet, Clade> getCoalescentCladeMap() {
//        return coalescentCladeMap;
//    }
//    
//    public Map<BitSet, Clade> getReticulationCladeMap() {
//        return reticulationCladeMap;
//    }

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
    

//    public double getSumCladeCredibility(NetworkNode node, BitSet bits) {
//
//        double sum = 0.0;
//
//        if (node.isLeaf()) {
//
//            int index = getTaxonIndex(node);
//            bits.set(2*index);
//        } else {
//
//            BitSet bits2 = new BitSet();
//            for (int i = 0; i < node.getChildCount(); i++) {
//
//            	NetworkNode node1 = node.getChild(i);
//
//                sum += getSumCladeCredibility(node1, bits2);
//            }
//
//            for (int i=1; i<bits2.length(); i=i+2) {
//                bits2.set(i, false);
//            }
//
//            if (node.isFake() && processSA) {
//                int index = getTaxonIndex(node.getDirectAncestorChild());
//                bits2.set(2 * index + 1);
//            }
//
//            sum += getCladeCredibility(bits2);
//
//            if (bits != null) {
//                bits.or(bits2);
//            }
//        }
//
//        return sum;
//    }
//
    
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

//    public BitSet removeClades(NetworkNode node, boolean includeTips) {
//
//        BitSet bits = new BitSet();
//
//        if (node.isLeaf()) {
//
//            int index = getTaxonIndex(node);
//            bits.set(2*index);
//
//            if (includeTips) {
//                removeClade(bits);
//            }
//
//        } else {
//
//            for (int i = 0; i < node.getChildCount(); i++) {
//
//            	NetworkNode node1 = node.getChild(i);
//
//                bits.or(removeClades(node1, includeTips));
//            }
//
//            for (int i=1; i<bits.length(); i=i+2) {
//                bits.set(i, false);
//            }
//            if (node.isFake() && processSA) {
//                int index = getTaxonIndex(node.getDirectAncestorChild());
//                bits.set(2 * index + 1);
//            }
//
//            removeClade(bits);
//        }
//
//        return bits;
//    }
//
//    private void removeClade(BitSet bits) {
//        Clade clade = cladeMap.get(bits);
//        if (clade != null) {
//            clade.setCount(clade.getCount() - 1);
//        }
//
//    }

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
    
    /**
     * clade needed to temporarily store reassortment clades before putting them all together
     *
     */
    public class DummyClade {
    	public DummyClade(BitSet bits, double height, boolean isLeft) {
	        this.bits = bits;
	        this.height = height;
	        this.isLeft = isLeft;
    	}
        
        public List<Object[]> getAttributeValues() {
            return attributeValues;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final DummyClade clade = (DummyClade) o;
            
            return !(bits != null ? !bits.equals(clade.bits) : clade.bits != null);

        }

        @Override
        public int hashCode() {
            return (bits != null ? bits.hashCode() : 0);
        }

        @Override
        public String toString() {
            return "clade1 " + bits + " height " + height;
//            return "count " + getCount();
       }
        
        BitSet bits;
        double height;
        boolean isLeft;
        List<Object[]> attributeValues = null;
    }


    public class ReassortmentClade {
        public ReassortmentClade(BitSet[] bits) {
            this.bits = Arrays.copyOf(bits, bits.length);
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

            final ReassortmentClade clade = (ReassortmentClade) o;
            

            for (int i = 0; i < bits.length; i++)
            	if (!(bits[i] == null && clade.bits[i] == null) || !(bits[i].equals(clade.bits[i])) )
            		return false;

            
            return true;

        }

        @Override
        public int hashCode() {
            return (bits != null ? bits.hashCode() : 0);
        }

        @Override
        public String toString() {
//            return "clade1 " + bits + " relative" + Arrays.toString(split) + " #" + count + " count " + getCount();
            return "count " + getCount();
       }

        int count;
        double credibility;
        BitSet[] bits;
        int[] split;
        List<Object[]> attributeValues = null;
    }
    
    
    public class BitSetArray {
        public BitSetArray(BitSet[] bits) {
            this.bits = Arrays.copyOf(bits, bits.length);
        }


        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final BitSetArray array = (BitSetArray) o;
        	return Arrays.equals(this.bits, array.bits);
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(bits);
        }

        BitSet[] bits;
    }

    
}