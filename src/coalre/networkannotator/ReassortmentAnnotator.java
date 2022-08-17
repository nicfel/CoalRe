/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package coalre.networkannotator;

import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A rewrite of TreeAnnotator that outputs how often reassortment events happen on trunk branches vs. other branches 
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class ReassortmentAnnotator {

	/**
	 * performs all the removing things steps
	 * @param network
	 * @param segmentToRemove
	 */
	void pruneNetwork(Network network, int[] segmentToRemove){
    	// remove all parts of the network that aren't informed by the genetic data
    	removeNonGeneticSegmentEdges(network);
    	for (int i = 0; i < segmentToRemove.length; i++)
    		removeSegment(network, segmentToRemove[i]);

    	// remove all loops
    	removeLoops(network);
    	// remove all empty edges in the segment
    	removeEmptyNetworkEdge(network);      
    	// remove all loops
    	removeLoops(network);

    	// remove all empty edges in the segment
    	removeEmptyNetworkEdge(network);  
	}
	
    /**
     * removes segments from network edges for which there is no
     * genetic information, i.e. segment that are above the segment tree root.
     * @param network
     */
    private void removeNonGeneticSegmentEdges(Network network){
    	// remove segments from edges if they are "above" the segment tree root
    	for (int i = 0; i < network.getSegmentCount(); i++){
    		removeSegmentFromEdge(network.getRootEdge(), i);
    	}    	
    }
    
    /**
     * remove segment from edges and child edges until a coalescent event with both children
     * carrying the segment is reached
     * @param edge
     * @param segIdx
     */
    private void removeSegmentFromEdge(NetworkEdge edge, int segIdx){
    	// remove the segment from the edge
    	edge.hasSegments.set(segIdx, false);
    	
    	// get all child segments
    	List<NetworkEdge> childEdges = edge.childNode.getChildEdges();
    	if (childEdges.size()==1){
    		removeSegmentFromEdge(childEdges.get(0), segIdx);
    	}else if (childEdges.size()==2){
    		// check if both children carry the segment
    		int carriesSegment = 0;
    		for (NetworkEdge childEdge : childEdges)
    			if (childEdge.hasSegments.get(segIdx))
    				carriesSegment++;
    		
    		// if it carries the segment, return and do nothing
    		if (carriesSegment==2){
    			return;
    		}else if (carriesSegment==1){
        		for (NetworkEdge childEdge : childEdges)
        			if (childEdge.hasSegments.get(segIdx))
        				removeSegmentFromEdge(childEdge, segIdx);

    		}else{
        		throw new IllegalArgumentException("at least one of the childre should carry the segment");  			
    		}		
    	}else{
    		throw new IllegalArgumentException("odd number of child edges");
    	}
    		
    }

    /**
     * removes all reticulation edges that start and end at the same place
     * @param network
     */
    private void removeLoops(Network network){
    	List<NetworkNode> reticulationNodes = network.getNodes().stream()
                .filter(e -> e.isReassortment())
                .filter(e -> e.getParentEdges().get(0).parentNode.equals(e.getParentEdges().get(1).parentNode))
                .filter(e -> e.getParentEdges().get(1).hasSegments.cardinality()>0)
                .collect(Collectors.toList());
    	
    	// for each of these, check if the parents are the same node
    	while (!reticulationNodes.isEmpty()){
    		NetworkNode node = reticulationNodes.get(0);
			// if this is the case, put all segment from 1 onto 0
			for (int i = 0; i < network.getSegmentCount(); i++){
				if (node.getParentEdges().get(1).hasSegments.get(i)){
					node.getParentEdges().get(0).hasSegments.set(i, true);
					node.getParentEdges().get(1).hasSegments.set(i, false);
				}    					
			}
			
			removeEmptyNetworkEdge(network);
			reticulationNodes = network.getNodes().stream()
	                .filter(e -> e.isReassortment())
	                .filter(e -> e.getParentEdges().get(0).parentNode.equals(e.getParentEdges().get(1).parentNode))
	                .filter(e -> e.getParentEdges().get(1).hasSegments.cardinality()>0)
	                .collect(Collectors.toList());			
    	}
    }

    /**
     * removes segment with id segIdx from the network.
     * @param network
     * @param segIdx
     */
    private void removeSegment(Network network, int segIdx){
    	// get all networkNodes
    	Set<NetworkEdge> networkEdges  = network.getEdges();
    	
    	// set carries segment nr segIdx to false for every node
    	for (NetworkEdge edge : networkEdges){
    		edge.hasSegments.set(segIdx, false);
    	}
    }
    
    /**
     * removes all edges from the network that don't carry any segments
     * @param network
     */
    private void removeEmptyNetworkEdge(Network network){
        List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        List<NetworkEdge> removableEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.hasSegments.cardinality()==0)
                .filter(e -> e.parentNode.isCoalescence())
                .collect(Collectors.toList());
        
        while (removableEdges.size()>0){
            int edgeInd = Randomizer.nextInt(removableEdges.size());     
            
        	removeEmptyReassortmentEdge(network, removableEdges.get(edgeInd));
        	
        	networkEdges = new ArrayList<>(network.getEdges());
            
            removableEdges = networkEdges.stream()
                    .filter(e -> !e.isRootEdge())
                    .filter(e -> e.childNode.isReassortment())
                    .filter(e -> e.hasSegments.cardinality()==0)
                    .filter(e -> e.parentNode.isCoalescence())
                    .collect(Collectors.toList());            
        } 

    }    

    private void removeEmptyReassortmentEdge(Network network, NetworkEdge edgeToRemove) {

        NetworkNode nodeToRemove = edgeToRemove.childNode;
        NetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
        NetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;

        // Remove edge and associated nodes
        NetworkEdge edgeToExtend = nodeToRemove.getChildEdges().get(0);
        nodeToRemove.removeChildEdge(edgeToExtend);
        nodeToRemove.removeParentEdge(edgeToRemove);
        nodeToRemove.removeParentEdge(edgeToRemoveSpouse);
        edgeToRemoveSpouseParent.removeChildEdge(edgeToRemoveSpouse);
        edgeToRemoveSpouseParent.addChildEdge(edgeToExtend);

        NetworkNode secondNodeToRemove = edgeToRemove.parentNode;
        NetworkEdge secondEdgeToExtend = getSisterEdge(edgeToRemove);

        secondNodeToRemove.removeChildEdge(secondEdgeToExtend);
        secondNodeToRemove.removeChildEdge(edgeToRemove);

        if (secondNodeToRemove.getParentEdges().get(0).isRootEdge()) {
            network.setRootEdge(secondEdgeToExtend);

        } else {
            NetworkEdge secondNodeToRemoveParentEdge = secondNodeToRemove.getParentEdges().get(0);
            NetworkNode secondNodeToRemoveParent = secondNodeToRemoveParentEdge.parentNode;
            secondNodeToRemoveParent.removeChildEdge(secondNodeToRemoveParentEdge);
            secondNodeToRemove.removeParentEdge(secondNodeToRemoveParentEdge);

            secondNodeToRemoveParent.addChildEdge(secondEdgeToExtend);
        }
    } 

    /**
     * Retrieve sister of given edge
     * @param childEdge child edge
     * @return sister of given child edge
     */
    private NetworkEdge getSisterEdge(NetworkEdge childEdge) {
        int idx = childEdge.parentNode.getChildEdges().indexOf(childEdge);
        int otherIdx = (idx + 1) % 2;

        return childEdge.parentNode.getChildEdges().get(otherIdx);
    }

    /**
     * Retrieve spouse of given edge
     * @param parentEdge parent edge
     * @return spouse of given parent edge
     */
    private NetworkEdge getSpouseEdge(NetworkEdge parentEdge) {
        int idx = parentEdge.childNode.getParentEdges().indexOf(parentEdge);
        int otherIdx = (idx + 1) % 2;

        return parentEdge.childNode.getParentEdges().get(otherIdx);
    }
    
}