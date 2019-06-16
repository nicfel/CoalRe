package coalre.operators;

import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.BitSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class DivertSegmentOperator extends EmptyEdgesNetworkOperator {

    @Override
    public double networkProposal() {
        double logHR = 0.0;

//        System.out.println(network.getExtendedNewick());

        List<NetworkEdge> sourceEdges = network.getEdges().stream()
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.hasSegments.cardinality()>0)
                .collect(Collectors.toList());

        if (sourceEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        logHR -= Math.log(1.0/sourceEdges.size());

        NetworkEdge sourceEdge = sourceEdges.get(Randomizer.nextInt(sourceEdges.size()));
        NetworkEdge destEdge = getSpouseEdge(sourceEdge);
        
//        if (sourceEdge.hasSegments.cardinality()==0)
//        	return Double.NEGATIVE_INFINITY;
        
        BitSet segsToDivert = getRandomUnconditionedSubset(sourceEdge.hasSegments);
        logHR -= getLogUnconditionedSubsetProb(sourceEdge.hasSegments);

//        System.out.println(sourceEdge.hasSegments);
//        System.out.println(segsToDivert);
//        System.exit(0);
        
        if (segsToDivert.cardinality()==0)
        	return Double.NEGATIVE_INFINITY;

        network.startEditing(this);
        

        logHR -= addSegmentsToAncestors(destEdge, segsToDivert);
        logHR += removeSegmentsFromAncestors(sourceEdge, segsToDivert);

        logHR += getLogUnconditionedSubsetProb(destEdge.hasSegments);

        int reverseSourceEdgeCount = (int)(network.getEdges().stream()
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.hasSegments.cardinality()>0)
                .count());

        logHR += Math.log(1.0/reverseSourceEdgeCount);

//        System.out.println(network.getExtendedNewick());
//        System.exit(0);
        return logHR;
    }


    /**
     * Remove segments from this edge and ancestors.
     *
     * @param edge edge at which to start removal
     * @param segsToRemove segments to remove from edge and ancestors
     * @return log probability of reverse operation
     */
    double removeSegmentsFromAncestors(NetworkEdge edge, BitSet segsToRemove) {
        double logP = 0.0;

        segsToRemove = (BitSet)segsToRemove.clone();
        segsToRemove.and(edge.hasSegments);

        if (segsToRemove.isEmpty())
            return logP;

        edge.hasSegments.andNot(segsToRemove);

        if (edge.isRootEdge())
            return logP;

        if (edge.parentNode.isReassortment()) {

            logP += Math.log(0.5)*segsToRemove.cardinality();

            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove);
            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(1), segsToRemove);

        } else {

            segsToRemove.andNot(getSisterEdge(edge).hasSegments);
            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove);

        }

        return logP;
    }

    /**
     * Add segments to this edge and ancestors.
     *
     * @param edge edge at which to start addition
     * @param segsToAdd segments to add to the edge and ancestors
     * @return log probability of operation
     */
    double addSegmentsToAncestors(NetworkEdge edge, BitSet segsToAdd) {
        double logP = 0.0;

        segsToAdd = (BitSet)segsToAdd.clone();
        segsToAdd.andNot(edge.hasSegments);

        if (segsToAdd.isEmpty())
            return logP;

        edge.hasSegments.or(segsToAdd);

        if (edge.isRootEdge())
            return logP;

        if (edge.parentNode.isReassortment()) {

            BitSet segsToAddLeft = new BitSet();
            BitSet segsToAddRight = new BitSet();

            for (int segIdx=segsToAdd.nextSetBit(0); segIdx != -1;
                    segIdx=segsToAdd.nextSetBit(segIdx+1)) {
                if (Randomizer.nextBoolean())
                    segsToAddLeft.set(segIdx);
                else
                    segsToAddRight.set(segIdx);

                logP += Math.log(0.5);
            }

            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAddLeft);
            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(1), segsToAddRight);

        } else {

            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAdd);
        }

        return logP;
    }
    
    protected BitSet getRandomUnconditionedSubset(BitSet sourceSegments) {
        BitSet destSegments = new BitSet();

        destSegments.clear();

        for (int segIdx = sourceSegments.nextSetBit(0); segIdx != -1;
             segIdx = sourceSegments.nextSetBit(segIdx + 1)) {

            if (Randomizer.nextBoolean())
                destSegments.set(segIdx);
        }

        return destSegments;
    }

    protected double getLogUnconditionedSubsetProb(BitSet sourceSegments) {
        return sourceSegments.cardinality()*Math.log(0.5);
    }


}
