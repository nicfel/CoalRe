package coalre.distribution.marginal;

import beast.core.Function;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.*;

class ParticleState {
    List<NetworkEdge> lineages;
    List<Node[]> segmentNodeArrays;

    int nSegments;

    public ParticleState(int nSegments) {
        this.nSegments = nSegments;
    }

    void clear() {
        lineages.clear();
        segmentNodeArrays.clear();
    }

    int size() {
        return lineages.size();
    }

    int getLineageCount() {
        return lineages.size();
    }

    NetworkEdge getNetworkLineage(int idx) {
        return lineages.get(idx);
    }

    Node[] getSegmentTreeNodes(int idx) {
        return segmentNodeArrays.get(idx);
    }

    void addLineage(NetworkEdge edge, Node[] segTreeNodeArray) {
        lineages.add(edge);
        segmentNodeArrays.add(segTreeNodeArray);
    }

    /**
     * Produce a deep copy of a given particle state.
     *
     * @param sourceParticle particle to copy
     * @param destParticle list to populate with copy
     * @return deep copy of state
     */
    void copyParticleState(ParticleState destParticle) {

        destParticle.clear();

        Map<NetworkNode,NetworkNode> nodeMap = new HashMap<>();

        for (int idx=0; idx<getLineageCount(); idx++) {
            destParticle.addLineage(getNetworkLineage(idx).getCopy(nodeMap),
                    Arrays.copyOf(getSegmentTreeNodes(idx),
                            getSegmentTreeNodes(idx).length));
        }
    }

    /**
     * Propagate network betwen observations
     *
     * @param particleState
     * @return
     */
    double propagateParticleState(double startTime, double nextObservationTime,
                                  Function reassortmentRate, Function populationSize) {
        double logConditionalP = 0.0;

        double currentTime = startTime;
        while (true) {
            int k = size();

            List<LineageIdxPair> coalescibleLineagePairs = getCoalescibleLineagePairs();
            int allowed_coalescences = coalescibleLineagePairs.size();
            int forbidden_coalescences = size()*(size()-1)/2 - allowed_coalescences;

            double a_coal_allowed = allowed_coalescences/populationSize.getArrayValue();
            double a_coal_forbidden = forbidden_coalescences/populationSize.getArrayValue();
            double a_reass = k*reassortmentRate.getArrayValue();

            double a_tot_allowed = a_coal_allowed + a_reass;
            double a_tot_forbidden = a_coal_forbidden;

            double nextEventTime = currentTime + Randomizer.nextExponential(a_tot_allowed);

            // Update particle weight
            double deltaT = Math.min(nextEventTime, nextObservationTime) - currentTime;
            logConditionalP += -deltaT * a_tot_forbidden;

            if (nextEventTime>nextObservationTime)
                break;

            // Implement event

            if (Randomizer.nextDouble()*a_tot_allowed < a_coal_allowed)
                coalesce(nextEventTime, coalescibleLineagePairs);
            else
                reassort(nextEventTime);

            currentTime = nextEventTime;
        }

        return logConditionalP;
    }

    private class LineageIdxPair {
        int idx1, idx2;

        public LineageIdxPair(int idx1, int idx2) {
            this.idx1 = idx1;
            this.idx2 = idx2;
        }
    }

   List<LineageIdxPair> getCoalescibleLineagePairs() {
       List<LineageIdxPair> idxPairs = new ArrayList<>();

       int allowed_coalescences = 0, forbidden_coalescences = 0;
       for (int i = 0; i < size(); i++) {
           NetworkEdge lineageA = getNetworkLineage(i);

           for (int j = 0; j < i; j++) {
               NetworkEdge lineageB = getNetworkLineage(j);

               if (lineageA.hasSegments.intersects(lineageB.hasSegments)) {
                   forbidden_coalescences += 1;
               } else {
                   allowed_coalescences += 1;
                   idxPairs.add(new LineageIdxPair(i, j));
               }
           }
       }

       return idxPairs;
   }

    void reassort(double time) {
        int idx = Randomizer.nextInt(getLineageCount());
        NetworkEdge lineage = getNetworkLineage(idx);
        Node[] segmentNodes = getSegmentTreeNodes(idx);

        BitSet hasSegs_left = new BitSet();
        BitSet hasSegs_right = new BitSet();

        for (int segIdx = lineage.hasSegments.nextSetBit(0);
             segIdx != -1; segIdx = lineage.hasSegments.nextSetBit(segIdx+1)) {
            if (Randomizer.nextBoolean()) {
                hasSegs_left.set(segIdx);
            } else {
                hasSegs_right.set(segIdx);
            }
        }

        // Stop here if reassortment event is unobservable
        if (hasSegs_left.cardinality() == 0 || hasSegs_right.cardinality() == 0)
            return;

        // Create reassortment node
        NetworkNode networkNode = new NetworkNode();
        networkNode.setHeight(time).addChildEdge(lineage);

        // Create reassortment lineages
        NetworkEdge leftLineage = new NetworkEdge(null, networkNode, hasSegs_left);
        Node[] leftSegNodes = new Node[nSegments];
        for (int segIdx = hasSegs_left.nextSetBit(0);
             segIdx != -1; segIdx = hasSegs_left.nextSetBit(segIdx+1)) {
            leftSegNodes[segIdx] = segmentNodes[segIdx];
        }

        NetworkEdge rightLineage = new NetworkEdge(null, networkNode, hasSegs_right);
        Node[] rightSegNodes = new Node[nSegments];
        for (int segIdx = hasSegs_right.nextSetBit(0);
             segIdx != -1; segIdx = hasSegs_right.nextSetBit(segIdx+1)) {
            rightSegNodes[segIdx] = segmentNodes[segIdx];
        }

        networkNode.addParentEdge(leftLineage);
        networkNode.addParentEdge(rightLineage);

        lineages.remove(idx);
        lineages.add(leftLineage);
        lineages.add(rightLineage);

        segmentNodeArrays.remove(idx);
        segmentNodeArrays.add(leftSegNodes);
        segmentNodeArrays.add(rightSegNodes);
    }

    private void coalesce(double coalescentTime, List<LineageIdxPair> coalescibleLineagePairs) {
        // Select lineages to coalesce

        LineageIdxPair coalescingPair = coalescibleLineagePairs.get(
                Randomizer.nextInt(coalescibleLineagePairs.size()));

        NetworkEdge lineageLeft = getNetworkLineage(coalescingPair.idx1);
        Node[] segNodesLeft = getSegmentTreeNodes(coalescingPair.idx1);
        NetworkEdge lineageRight = getNetworkLineage(coalescingPair.idx2);
        Node[] segNodesRight = getSegmentTreeNodes(coalescingPair.idx2);

        // Create coalescent node
        NetworkNode coalescentNode = new NetworkNode();
        coalescentNode.setHeight(coalescentTime)
                .addChildEdge(lineageLeft)
                .addChildEdge(lineageRight);
        lineageLeft.parentNode = coalescentNode;
        lineageRight.parentNode = coalescentNode;

        // Merge segment flags:
        BitSet hasSegments = new BitSet();
        hasSegments.or(lineageLeft.hasSegments);
        hasSegments.or(lineageRight.hasSegments);

        // Merge segment node arrays:

        Node[] segmentNodes = new Node[nSegments];
        for (int segIdx = lineageLeft.hasSegments.nextSetBit(0);
             segIdx != -1; segIdx = lineageLeft.hasSegments.nextSetBit(segIdx+1)) {
            segmentNodes[segIdx] = segNodesLeft[segIdx];
        }

        for (int segIdx = lineageRight.hasSegments.nextSetBit(0);
             segIdx != -1; segIdx = lineageRight.hasSegments.nextSetBit(segIdx+1)) {
            segmentNodes[segIdx] = segNodesRight[segIdx];
        }

        // Create new lineage
        NetworkEdge lineage = new NetworkEdge(null, coalescentNode, hasSegments);
        coalescentNode.addParentEdge(lineage);

        lineages.remove(lineageLeft);
        lineages.remove(lineageRight);
        lineages.add(lineage);

        segmentNodeArrays.remove(segNodesLeft);
        segmentNodeArrays.remove(segNodesRight);
        segmentNodeArrays.add(segmentNodes);
    }

    double getObservedEventProbability(ObservedEvent event,
                                       Function reassortmentRate,
                                       Function populationSize) {

        switch(event.type) {
            case SAMPLE:
            case COALESCENCE:
        }

        return 0.0;
    }

}
