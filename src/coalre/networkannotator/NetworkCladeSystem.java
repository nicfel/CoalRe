package coalre.networkannotator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import javax.swing.plaf.synth.SynthSeparatorUI;

import beast.base.util.DiscreteStatistics;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;


/**
 * extracted from TreeAnnotator
 */
public class NetworkCladeSystem {

	
    protected List<Map<Double, DummyClade>> newCoalescentCladeMap = new ArrayList<>();
    protected Map<BitSetArray, ReassortmentClade> cladeMap = new HashMap<>();
    protected List<Map<Double, DummyClade>> newReassortmentCladeMap = new ArrayList<>();
    protected Map<BitSetArray, ReassortmentClade> reassortmentCladeMap = new HashMap<>();
    public List<String> leafNodeMap;
    
    protected Map<BitSet, DummyClade> leafCladeMap = new HashMap<>();


    public boolean[] followSegmentAlready;
    protected boolean started;
    public NetworkCladeSystem() { 
    }
    
    /**
     * adds all leaf labels to a list
     * @param leafNodes
     * @param nrSegments
     */
    public void setLeafLabels(List<NetworkNode> leafNodes, int nrSegments){
    	cladeMap = new HashMap<>();
    	reassortmentCladeMap = new HashMap<>();
    	leafCladeMap = new HashMap<>();
    	
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
    	newCoalescentCladeMap = new ArrayList<>();
		newReassortmentCladeMap = new ArrayList<>();
		for (int i = 0; i < network.getSegmentCount(); i++){
			newCoalescentCladeMap.add(new HashMap<>());
			newReassortmentCladeMap.add(new HashMap<>());			
			final int segment = i;
			// get the segment root node
	        List<NetworkNode> rootEdge = network.getNodes().stream()
	                .filter(e -> e.isCoalescence())
	                .filter(e -> !e.getParentEdges().get(0).hasSegments.get(segment))
	                .filter(e -> e.getChildEdges().get(0).hasSegments.get(segment))
	                .filter(e -> e.getChildEdges().get(1).hasSegments.get(segment))
	                .collect(Collectors.toList());
	        
	        
	        if (rootEdge.size()==1)
	        	addClades(rootEdge.get(0), includeTips, i);
	        else if (rootEdge.size()>1)
	        	throw new IllegalArgumentException("segment tree root not found");	        
		}
		
		// build the reassortment clades with all segments
		buildCoalescentCladeMap(network.getSegmentCount());		
		buildReassortmentCladeMap(network.getSegmentCount());	
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
            for (NetworkEdge childEdge : childEdges){
	        	if (childEdge.hasSegments.get(segment))
	    			bits.or(addClades(childEdge.childNode, includeTips, segment));
            	
            }
            
            // if node is coalescent, add the bitset if the coalescent event is observed on a segment tree
            if (node.isCoalescence()){
        		if (node.getChildEdges().get(0).hasSegments.get(segment) &&
        			node.getChildEdges().get(1).hasSegments.get(segment)){
        			addCoalescentClade(bits, segment, node.getHeight());
        		}    			
            }else{
            	if (node.getParentEdges().get(0).hasSegments.get(segment))
            		addReassortmentClade(bits, segment, node.getHeight(), false);
            	else if (node.getParentEdges().get(1).hasSegments.get(segment))
            		addReassortmentClade(bits, segment, node.getHeight(), true);
            }     
            
        }

        return bits;
    }

    /**
     * adds a new dummy clade to a coalescent node in the network
     * @param bits
     * @param segment
     * @param height
     */
    private void addCoalescentClade(BitSet bits, int segment, Double height) {
    	DummyClade clade = newCoalescentCladeMap.get(segment).get(height);
        if (clade == null) {
            clade = new DummyClade(bits, height, false);
            newCoalescentCladeMap.get(segment).put(height, clade);
        }else{
        	throw new IllegalArgumentException("coalescent clade should never be found");
        }
    }    
    
    /**
     * adds new dummy clades to the new clade map
     * @param bits
     * @param segment
     * @param height
     */
    protected void addReassortmentClade(BitSet bits, int segment, Double height, boolean isLeft) {
        DummyClade clade = newReassortmentCladeMap.get(segment).get(height);
        if (clade == null) {
            clade = new DummyClade(bits, height, isLeft);
            newReassortmentCladeMap.get(segment).put(height, clade);
        }else{
        	throw new IllegalArgumentException("reassortment clade should never be found");
        }
    }   
    
    private void buildCoalescentCladeMap(int nrSegments) {
    	// get all unique reassortment node heights
    	List<Double> nodeHeights = new ArrayList<>();
    	for (int i = 0; i < newCoalescentCladeMap.size(); i++){
    		for (Double height : newCoalescentCladeMap.get(i).keySet()){
    			int index = nodeHeights.indexOf(height);
    			if (index==-1){
    				nodeHeights.add(height);    			
    			}
    		}
    	}    	
    	
    	for (int i = 0; i < nodeHeights.size(); i++){    	
    		// make a bit set array that is empty if a segment goes left
    		BitSet[] bits = new BitSet[nrSegments];
    		for (int j = 0; j < nrSegments; j++){
    			if (newCoalescentCladeMap.get(j).get(nodeHeights.get(i))!=null)
    				bits[j] = newCoalescentCladeMap.get(j).get(nodeHeights.get(i)).bits;    			
    		}
    		
    		BitSetArray bitsArray = new BitSetArray(bits);
    		
    		// add the bits to a new reassortmentClade
            ReassortmentClade clade = cladeMap.get(bitsArray);
            if (clade == null) {
                clade = new ReassortmentClade(bits);
            	cladeMap.put(bitsArray, clade);                
            }              
            clade.setCount(clade.getCount() + 1);
    	}    	
    }       
   
    protected void buildReassortmentCladeMap(int nrSegments) {
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

    public void collectAttributes(Network network, Set<String> attributeNames, boolean includeTips) {
		newCoalescentCladeMap = new ArrayList<>();
		newReassortmentCladeMap = new ArrayList<>();
		for (int i = 0; i < network.getSegmentCount(); i++){
			newCoalescentCladeMap.add(new HashMap<>());
			newReassortmentCladeMap.add(new HashMap<>());
			
			
			final int segment = i;
			// get the segment root node
	        List<NetworkNode> rootEdge = network.getNodes().stream()
	                .filter(e -> e.isCoalescence())
	                .filter(e -> !e.getParentEdges().get(0).hasSegments.get(segment))
	                .filter(e -> e.getChildEdges().get(0).hasSegments.get(segment))
	                .filter(e -> e.getChildEdges().get(1).hasSegments.get(segment))
	                .collect(Collectors.toList());
	        
	        
	        if (rootEdge.size()==1)
	        	collectAttributes(rootEdge.get(0), attributeNames, includeTips, i);
	        else if (rootEdge.size()>1)
	        	throw new IllegalArgumentException("segment tree root not found");	        
		}		

		collectAtributesCoalescentCladeMap(network.getSegmentCount(), attributeNames);
		collectAtributesReassortmentCladeMap(network.getSegmentCount(), attributeNames);		
    }
    
    private BitSet collectAttributes(NetworkNode node, Set<String> attributeNames, boolean includeTips, int segment) {
        BitSet bits = new BitSet();
        
        
        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            if (index < 0) {
                throw new IllegalArgumentException("Taxon with height= " + node.getHeight() + ", not found in target tree");
            }
        	if (node.getParentEdges().get(0).hasSegments.get(segment))
        			bits.set(index);
        	
        	if (includeTips && node.getParentEdges().get(0).hasSegments.get(segment))
    			collectAttributesForLeaves(bits, node, attributeNames, segment, node.getHeight());


        } else {        	
           	
        	// get the children of that node
            List<NetworkEdge> childEdges = node.getChildEdges();
            // add all children to the bitset
            for (NetworkEdge childEdge : childEdges){
	        	if (childEdge.hasSegments.get(segment))
	    			bits.or(collectAttributes(childEdge.childNode, attributeNames, includeTips, segment));
            }  
            
            if (node.isCoalescence()){
        		if(node.getChildEdges().get(0).hasSegments.get(segment) &&
        			node.getChildEdges().get(1).hasSegments.get(segment)){
        			collectAttributesForClade(bits, node, attributeNames, segment, node.getHeight());
        		}
            } else {
            	if (node.getParentEdges().get(0).hasSegments.get(segment))
            		collectAttributesForReassortmentClade(bits, node, attributeNames, segment, node.getHeight(), false);
            	else if (node.getParentEdges().get(1).hasSegments.get(segment))
            		collectAttributesForReassortmentClade(bits, node, attributeNames, segment, node.getHeight(), true);
            }
            
        }

        return bits;
    }
    
    private void collectAttributesForClade(BitSet bits, NetworkNode node, Set<String> attributeNames, int segment, Double height) {
        DummyClade clade = newCoalescentCladeMap.get(segment).get(height);
        if (clade == null) {
        	
        	clade = new DummyClade(bits, height, false);

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
//                        value = getBranchLength(node);
//                        break;
                    default:
                    	throw new IllegalArgumentException("Summary not implemented for values other than height and posterior");
                }

                values[i] = value;

                i++;
            }
            clade.attributeValues.add(values);
            
            newCoalescentCladeMap.get(segment).put(height, clade);

        }
    }
    
    private void collectAttributesForLeaves(BitSet bits, NetworkNode node, Set<String> attributeNames, int segment, Double height) {
    	DummyClade clade = leafCladeMap.get(bits);
        if (clade == null) {
        	
        	clade = new DummyClade(bits, height, false);

            clade.attributeValues = new ArrayList<>();
            
            int i = 0;
            Object[] values = new Object[attributeNames.size()];
            for (String attributeName : attributeNames) {

                Object value;
                switch (attributeName) {
                    case "height":
                        value = node.getHeight();
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
            
            leafCladeMap.put(bits, clade);
        }else{
            int i = 0;
            Object[] values = new Object[attributeNames.size()];
            for (String attributeName : attributeNames) {

                Object value;
                switch (attributeName) {
                    case "height":
                        value = node.getHeight();
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
//        	System.err.println("reassortment clade should never be found");
       	throw new IllegalArgumentException("reassortment clade should never be found");
        }
    }

    private void collectAtributesCoalescentCladeMap(int nrSegments, Set<String> attributeNames) {
    	// get all unique reassortment node heights
    	List<Double> nodeHeights = new ArrayList<>();
    	for (int i = 0; i < newCoalescentCladeMap.size(); i++){
    		for (Double height : newCoalescentCladeMap.get(i).keySet()){
    			int index = nodeHeights.indexOf(height);
    			if (index==-1){
    				nodeHeights.add(height);
    			}
    		}
    	}    
    	    	
    	for (int i = 0; i < nodeHeights.size(); i++){
    		// make a bit set array that is empty if a segment goes left
    		BitSet[] bits = new BitSet[nrSegments];
    		for (int j = 0; j < nrSegments; j++){
    			if (newCoalescentCladeMap.get(j).get(nodeHeights.get(i))!=null)
    				bits[j] = newCoalescentCladeMap.get(j).get(nodeHeights.get(i)).bits;    			
    		}
    		
    		
    		BitSetArray bitsArray = new BitSetArray(bits);
    		
    		// add the bits to a new reassortmentClade
            ReassortmentClade clade = cladeMap.get(bitsArray);
            if (clade != null) {
            	if (clade.attributeValues==null)
            		clade.attributeValues = new ArrayList<>();
            		
            	// add the attributes
            	for (int j = 0; j < nrSegments; j++){
            		if (newCoalescentCladeMap.get(j).get(nodeHeights.get(i))!=null){
            			clade.attributeValues.addAll(newCoalescentCladeMap.get(j).get(nodeHeights.get(i)).getAttributeValues());
            		}
            	}
            	clade.setCount(clade.getCount()+1);
            }
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
     
    public void summarizeAttributes(Network network, Set<String> attributeNames, boolean useMean, int nrNetworks, boolean onTarget) {
		boolean[] followSegment = new boolean[network.getSegmentCount()];
		for (int i=0;i<network.getSegmentCount();i++)followSegment[i] = false;
		
		followSegmentAlready = Arrays.copyOf(followSegment, followSegment.length);
		
		//
		Map<NetworkNode, BitSet[]> passedBitSet = new HashMap<>();
		
		
    	// summarizes all coalescent events
    	summarizeAttributes(network.getRootEdge().childNode, attributeNames, useMean, nrNetworks, network.getSegmentCount(), followSegment, onTarget, passedBitSet);    	
    }
    
    private BitSet[] summarizeAttributes(NetworkNode node, Set<String> attributeNames, boolean useMean, 
    		int nrNetworks, int nrSegments, boolean[] followSegment_in, boolean onTarget, Map<NetworkNode, BitSet[]> passedBitSet) {

		// check if the node was already passed
    	if (passedBitSet.get(node)!=null)
    		return passedBitSet.get(node);

    	
        BitSet[] bits = new BitSet[nrSegments];
        for (int i = 0; i < nrSegments; i++) bits[i] = new BitSet();
        
        boolean[] followSegment = Arrays.copyOf(followSegment_in, followSegment_in.length);

        if (node.isLeaf()) {
            int index = getTaxonIndex(node);
            if (index < 0) {
                throw new IllegalArgumentException("Taxon with height= " + node.getHeight() + ", not found in target tree");
            }
            for (int i = 0; i < nrSegments; i++){
            	if (node.getParentEdges().get(0).hasSegments.get(i))
            		bits[i].set(index);
            }
            summarizeAttributesForLeaf(bits, node, attributeNames, useMean, nrNetworks, nrSegments, onTarget);
        } else {
        	
            // check if the node is the root of a segment tree, if so, follow the segment
            if (node.isCoalescence()){
	           	 for (int i = 0; i < nrSegments; i++){
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
				
				BitSet[] newBits = summarizeAttributes(childEdge.childNode, attributeNames, useMean, nrNetworks, nrSegments, followSegmentout, onTarget, passedBitSet);
            	for (int i = 0; i < nrSegments;i++){
	            	if (childEdge.hasSegments.get(i)){
            			bits[i].or(newBits[i]);	   
	            	}
            	}
            	passedBitSet.put(node, bits);
				
            }            
            if (!node.isReassortment())
            	summarizeAttributesForClade(bits, node, attributeNames, useMean, nrNetworks, nrSegments, onTarget);   
            else{
            	summarizeAttributesForReassortmentClade(bits, node, attributeNames, useMean, nrNetworks, nrSegments, onTarget);
            }
        }  
        
        return bits;
    }
    
    private void summarizeAttributesForClade(BitSet[] bits, NetworkNode node, 
    		Set<String> attributeNames, boolean useMean, int nrNetworks, int nrSegments, boolean onTarget) {

		BitSet[] newBits = new BitSet[nrSegments];
		for (int i = 0 ; i < newBits.length; i++)
			if (bits[i].cardinality()>0 && 
					node.getChildEdges().get(0).hasSegments.get(i) &&
					node.getChildEdges().get(1).hasSegments.get(i)) newBits[i] = bits[i];
		
		BitSetArray bitsArray = new BitSetArray(newBits);
    	
    	// check if every bitset element is null, can happen if segments are removed
    	boolean allNull = true;
		for (int i = 0; i < bits.length; i++)
			if (newBits[i]!=null)
				allNull = false;
		
		// keeps track of the height of the target network,
		Double targetHeight = null;

    	
    	if (!allNull){	    	
	        for (String attributeName : attributeNames) {
	            switch (attributeName) {
	            	case "height":
	            		node.setMetaData("");
	            			            		
	            		List<Double> height = new ArrayList<>();
	            		if(cladeMap.get(bitsArray).count>1){
	            			// can happend when a target network is used
				    		List<Object[]> rawHeights = cladeMap.get(bitsArray).getAttributeValues();
				            for (int i = 0; i < rawHeights.size(); i++)
				            	height.add((double) rawHeights.get(i)[0]);
	            		}

			            
	            		if(onTarget)
	            			targetHeight = node.getHeight();            	
			            
	            		double posterior = (double) (cladeMap.get(bitsArray).getCount()-1)/ (double) (nrNetworks);
	            			            		
	            		// Convert height to Array
	            		double[] heightarray = new double[height.size()];
	            		for (int i = 0; i < height.size(); i++)
	            			heightarray[i] = height.get(i);
	            		
	            		if(!onTarget){
		            		if (heightarray.length>0){
								if (useMean){
									node.setHeight(DiscreteStatistics.mean(heightarray));
								}else{
								    node.setHeight(DiscreteStatistics.median(heightarray));
								}
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
	            		
	            		if (targetHeight!=null){
	            			node.setMetaData(",posterior=" + posterior +
	            					",targetHeight=" + targetHeight +
	            					",height_95%_HPD={" + minHPD + "," + maxHPD + "}" + 
			            		node.getMetaData() + "");

	            		}else{
	            			node.setMetaData(",posterior=" + posterior +
	            					",height_95%_HPD={" + minHPD + "," + maxHPD + "}" + 
			            		node.getMetaData() + "");
	            		}

	                   
	                case "length":
	                    break;
	                default:
	                	throw new IllegalArgumentException("");
	            }
	        }
    	}
    	
    }
    
    private void summarizeAttributesForLeaf(BitSet[] bits, NetworkNode node, 
    		Set<String> attributeNames, boolean useMean, int nrNetworks, int nrSegments, boolean onTarget) {

		BitSet newBits = new BitSet();
		for (int i = 0 ; i < bits.length; i++)
			if (bits[i].cardinality()>0 && 
					node.getParentEdges().get(0).hasSegments.get(i)) newBits = bits[i];
		
		
		// keeps track of the height of the target network,
		Double targetHeight = null;

    	
        for (String attributeName : attributeNames) {
            switch (attributeName) {
            	case "height":
            		node.setMetaData("");
            			            		
            		List<Double> height = new ArrayList<>();
		    		List<Object[]> rawHeights = leafCladeMap.get(newBits).getAttributeValues();
		            for (int i = 0; i < rawHeights.size(); i++)
		            	height.add((double) rawHeights.get(i)[0]);

		            
            		if(onTarget)
            			targetHeight = node.getHeight();            	
            			            		
            		// Convert height to Array
            		double[] heightarray = new double[height.size()];
            		for (int i = 0; i < height.size(); i++)
            			heightarray[i] = height.get(i);
            		
            		if(!onTarget){
	            		if (heightarray.length>0){
							if (useMean){
								node.setHeight(DiscreteStatistics.mean(heightarray));
							}else{
							    node.setHeight(DiscreteStatistics.median(heightarray));
							}
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
            		
            		if (targetHeight!=null){
            			node.setMetaData(",targetHeight=" + targetHeight +
            					",height_95%_HPD={" + minHPD + "," + maxHPD + "}" + 
		            		node.getMetaData() + "");

            		}else{
            			node.setMetaData(",height_95%_HPD={" + minHPD + "," + maxHPD + "}" + 
		            		node.getMetaData() + "");
            		}

                   
                case "length":
                    break;
                default:
                	throw new IllegalArgumentException("");
            }
        }
    	
    }
        
    private void summarizeAttributesForReassortmentClade(BitSet[] bits, NetworkNode node, 
    		Set<String> attributeNames, boolean useMean, int nrNetworks, int nrSegments, boolean onTarget) {
    	
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
    	
    	// check if every bitset element is null, can happen if segments are removed
    	boolean allNull = true;
		for (int i = 0; i < bitsarray.length; i++)
			if (bitsarray[i]!=null)
				allNull = false;
		
		// keeps track of the height of the target network,
		Double targetHeight = null;
    	
    	if (!allNull){
	        for (String attributeName : attributeNames) {
	            switch (attributeName) {
	            	case "height":
	            		node.setMetaData("");
	            		
	            		List<Double> height = new ArrayList<>();
	            		
	            		if(reassortmentCladeMap.get(keyArray).count>1){
	            			// can happend when a target network is used
				    		List<Object[]> rawHeights = reassortmentCladeMap.get(keyArray).getAttributeValues();
				            for (int i = 0; i < rawHeights.size(); i++)
				            	height.add((double) rawHeights.get(i)[0]);
	            		}else if(!onTarget){
	            			throw new IllegalArgumentException("reassortment clade not found");
	            		}
	
	            		
	            		double posterior = (double) (reassortmentCladeMap.get(keyArray).getCount()-1)/ (double) (nrNetworks);
	            		
	            		if(onTarget)
	            			targetHeight = node.getHeight();
	            		
	            		// Convert height to Array
	            		double[] heightarray = new double[height.size()];
	            		for (int i = 0; i < height.size(); i++)
	            			heightarray[i] = height.get(i);
	            			            		
	            		if(!onTarget){
		            		if (heightarray.length>0){
								if (useMean){
									node.setHeight(DiscreteStatistics.mean(heightarray));
								}else{
								    node.setHeight(DiscreteStatistics.median(heightarray));
								}
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
	            		
	            		if (targetHeight!=null){
	            			node.setMetaData(",posterior=" + posterior +
	            					",targetHeight=" + targetHeight +
	            					",height_95%_HPD={" + minHPD + "," + maxHPD + "}" + 
			            		node.getMetaData() + "");

	            		}else{
	            			node.setMetaData(",posterior=" + posterior +
	            					",height_95%_HPD={" + minHPD + "," + maxHPD + "}" + 
			            		node.getMetaData() + "");
	            		}
	
	                case "length":
	                    break;
	                default:
	                	throw new IllegalArgumentException("");
	            }
	        }
    	}
    	
    }
    
    public void calculateCladeCredibilities(int totalTreesUsed) {
        for (ReassortmentClade clade : cladeMap.values()) {	        	
            if (clade.getCount() > totalTreesUsed) {
            	for (int j = 0; j < leafNodeMap.size();j++)
            		if (clade.bits[0].get(j))
            			System.err.println(leafNodeMap.get(j));

                throw new AssertionError("clade.getCount=(" + clade.getCount() +
                        ") should be <= totalTreesUsed = (" + totalTreesUsed + ")");
            }
            clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);
        }
        
        for (ReassortmentClade clade : reassortmentCladeMap.values()) {	        	
            if (clade.getCount() > totalTreesUsed) {
            	for (int j = 0; j < leafNodeMap.size();j++)
            		if (clade.bits[0].get(j))
            			System.err.println(leafNodeMap.get(j));

                throw new AssertionError("clade.getCount=(" + clade.getCount() +
                        ") should be <= totalTreesUsed = (" + totalTreesUsed + ")");
            }
            clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);
        }

    	
    }
        
    public double getLogCladeCredibility(Network network){
    	newCoalescentCladeMap = new ArrayList<>();
		newReassortmentCladeMap = new ArrayList<>();
		for (int i = 0; i < network.getSegmentCount(); i++){
			newCoalescentCladeMap.add(new HashMap<>());
			newReassortmentCladeMap.add(new HashMap<>());
			final int segment = i;
			// get the segment root node
	        List<NetworkNode> rootEdge = network.getNodes().stream()
	                .filter(e -> e.isCoalescence())
	                .filter(e -> !e.getParentEdges().get(0).hasSegments.get(segment))
	                .filter(e -> e.getChildEdges().get(0).hasSegments.get(segment))
	                .filter(e -> e.getChildEdges().get(1).hasSegments.get(segment))
	                .collect(Collectors.toList());
	        
	        if (rootEdge.size()==1)
	        	addClades(rootEdge.get(0), false, i);
	        else if (rootEdge.size()>1)
	        	throw new IllegalArgumentException("segment tree root not found");	        
		}	
		
		return computeSumLogCladeCredibility(network.getSegmentCount());

    }
    
    public double computeSumLogCladeCredibility(int nrSegments) {
    	
        double logCladeCredibility = 0.0;

    	// get all unique reassortment node heights
    	List<Double> nodeHeights = new ArrayList<>();
    	for (int i = 0; i < newCoalescentCladeMap.size(); i++){
    		for (Double height : newCoalescentCladeMap.get(i).keySet()){
    			int index = nodeHeights.indexOf(height);
    			if (index==-1){
    				nodeHeights.add(height);    			
    			}
    		}
    	}    	
    	
    	for (int i = 0; i < nodeHeights.size(); i++){   		
    		// make a bit set array that is empty if a segment goes left
    		BitSet[] bits = new BitSet[nrSegments];
    		for (int j = 0; j < nrSegments; j++){
    			if (newCoalescentCladeMap.get(j).get(nodeHeights.get(i))!=null)
    				bits[j] = newCoalescentCladeMap.get(j).get(nodeHeights.get(i)).bits;    			
    		}
    		
    		BitSetArray bitsArray = new BitSetArray(bits);
    		
    		// add the bits to a new reassortmentClade
            ReassortmentClade clade = cladeMap.get(bitsArray);
            if (clade == null) {
                throw new IllegalArgumentException("coalescence clade not found");              
            } else {
            	logCladeCredibility += Math.log(clade.getCredibility())*clade.getNotNull();
            }   
    	} 
    	
    	
        return logCladeCredibility;
    }
    
    /**
     * get the index of a leaf node
     */ 
    protected int getTaxonIndex(NetworkNode leaf){
    	return leafNodeMap.indexOf(leaf.getTaxonLabel());
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
            notnull = 0;
            for (int i = 0; i < bits.length; i++)
            	if (bits[i]!=null)
            		notnull++;
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
        
        public int getNotNull(){
        	return notnull;
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
        int notnull;
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