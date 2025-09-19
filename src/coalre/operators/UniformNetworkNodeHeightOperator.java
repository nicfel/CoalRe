package coalre.operators;

import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.List;
import java.util.stream.Collectors;

public class UniformNetworkNodeHeightOperator extends NetworkOperator {

    @Override
    public double networkProposal() {
        List<NetworkNode> networkNodes = network.getNodes().stream()
                .filter(n -> !n.getParentEdges().get(0).isRootEdge())
                
                .collect(Collectors.toList());

        if (networkNodes.isEmpty())
            return Double.NEGATIVE_INFINITY;

        NetworkNode node = networkNodes.get(Randomizer.nextInt(networkNodes.size()));

        double maxHeight = Double.POSITIVE_INFINITY;
        for (NetworkEdge parentEdge : node.getParentEdges())
            if (parentEdge.parentNode.getHeight() < maxHeight)
                maxHeight = parentEdge.parentNode.getHeight();

        double minHeight = Double.NEGATIVE_INFINITY;
        for (NetworkEdge childEdge : node.getChildEdges())
            if (childEdge.childNode.getHeight() > minHeight)
                minHeight = childEdge.childNode.getHeight();

        network.startEditing(this);
        node.setHeight(minHeight + Randomizer.nextDouble()*(maxHeight-minHeight));
        
        if (node.isCoalescence()) {
	        if (node.segmentIndices!=null && segmentTreesInput.get().size()>0) {
	        	for (int i =0; i < node.segmentIndices.length; i++) {
	        		if (node.getChildEdges().get(0).hasSegments.get(i) && node.getChildEdges().get(1).hasSegments.get(i))
	        			segmentTreesInput.get().get(i).getNode(node.segmentIndices[i]).setHeight(node.getHeight());
				}
	        }
        }
        segmentsChanged.set(0, network.getSegmentCount(), false);

        
        return 0.0;
    }

}
