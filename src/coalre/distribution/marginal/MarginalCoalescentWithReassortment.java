package coalre.distribution.marginal;

import beast.core.Distribution;
import beast.core.Function;
import beast.core.Input;
import beast.core.State;
import beast.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.*;

public class MarginalCoalescentWithReassortment extends Distribution {

    public Input<Function> reassortmentRateInput = new Input<>(
	        "reassortmentRate",
            "reassortment rate (per lineage per unit time)",
            Input.Validate.REQUIRED);

	public Input<Function> populationSizeInput = new Input<>(
	        "populationSize",
            "(Constant) size of population.",
            Input.Validate.REQUIRED);

	public Input<ObservedEventList> eventListInput = new Input<>(
	        "observedEventList",
            "List of observed segment tree events.",
            Input.Validate.REQUIRED);

    public Input<Integer> nParticlesInput = new Input<>(
            "nParticles",
            "Number of particles to use in SMC calculation.",
            Input.Validate.REQUIRED);

    Function reassortmentRate;
    Function populationSize;

    ObservedEventList eventList;

    int nParticles;
    List<NetworkEdge>[] particleStates, particleStatesPrime;
    double[] particleWeights;


    public MarginalCoalescentWithReassortment() { }

    @Override
    public void initAndValidate() {
        eventList = eventListInput.get();
        nParticles = nParticlesInput.get();

        reassortmentRate = reassortmentRateInput.get();
        populationSize = populationSizeInput.get();

        particleStates = new List[nParticles];
        for (int p=0; p<nParticles; p++) {
            particleStates[p] = new ArrayList<>();
            particleStatesPrime[p] = new ArrayList<>();
            particleWeights[p] = 0.0;
        }
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        // Initialize particle states

        for (int p=0; p<nParticles; p++)
            particleStates[p].clear();

        // Loop over observations

        double prevEventTime = 0.0;
        for (ObservedEvent event : eventList.getEventList()) {

            // Propagate particles

            double maxWeight = Double.NEGATIVE_INFINITY;
            for (int p=0; p<nParticles; p++) {

                particleWeights[p] = propagateParticleState(particleStates[p], prevEventTime, event.time);

                if (particleWeights[p]>Double.NEGATIVE_INFINITY)
                    particleWeights[p] += getObservedEventProbability(particleStates[p], event);

                if (particleWeights[p] > maxWeight)
                    maxWeight = particleWeights[p];
            }

            // Compute average weight

            double sumScaledWeights = 0.0;
            for (int p=0; p<nParticles; p++) {
                particleWeights[p] = particleWeights[p] - maxWeight;
                sumScaledWeights += Math.exp(particleWeights[p]);
            }

            logP += Math.log(sumScaledWeights/nParticles) + maxWeight;

            // Resample particles

            double[] u = new double[nParticles];
            for (int p=0; p<nParticles; p++)
                u[p] = Randomizer.nextDouble()*sumScaledWeights;

            Arrays.sort(u);

            double accumulator = 0;
            int pPrime = 0;
            for (int p=0; p<nParticles; p++) {

                while (u[pPrime]-accumulator < Math.exp(particleWeights[p])) {
                    copyParticleState(particleStates[p], particleStatesPrime[pPrime]);
                    pPrime += 1;
                }

                accumulator += Math.exp(particleWeights[p]);
            }

            particleStates = particleStatesPrime;

        }


        return logP;
    }

    /**
     * Produce a deep copy of a given particle state.
     *
     * @param sourceParticle particle to copy
     * @param destParticle list to populate with copy
     * @return deep copy of state
     */
    void copyParticleState(List<NetworkEdge> sourceParticle, List<NetworkEdge> destParticle) {

        destParticle.clear();

        Map<NetworkNode,NetworkNode> nodeMap = new HashMap<>();

        for (NetworkEdge edge : sourceParticle) {
            destParticle.add(edge.getCopy(nodeMap));
        }
    }

    /**
     * Propagate network betwen observations
     *
     * @param particleState
     * @return
     */
    double propagateParticleState(List<NetworkEdge> particleState, double startTime, double nextObservationTime) {
        double logConditionalP = 0.0;

        double currentTime = startTime;
        while (true) {
            int k = particleState.size();

            List<LineagePair> coalescibleLineagePairs = new ArrayList<>();

            int allowed_coalescences = 0, forbidden_coalescences = 0;
            for (int i=0; i<particleState.size(); i++) {
                NetworkEdge lineageA = particleState.get(i);

                for (int j=0; j<i; j++) {
                    NetworkEdge lineageB = particleState.get(i);

                    if (lineageA.hasSegments.intersects(lineageB.hasSegments)) {
                        forbidden_coalescences += 1;
                    } else {
                        allowed_coalescences += 1;
                        coalescibleLineagePairs.add(new LineagePair(lineageA, lineageB));
                    }
                }
            }

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
                coalesce(nextEventTime, particleState, coalescibleLineagePairs);
            else
                reassort(nextEventTime, particleState);

            currentTime = nextEventTime;
        }

        return logConditionalP;
    }

    class LineagePair {
        NetworkEdge lineage1, lineage2;

        public LineagePair(NetworkEdge lineage1, NetworkEdge lineage2) {
            this.lineage1 = lineage1;
            this.lineage2 = lineage2;
        }
    }

    private void coalesce(double coalescentTime, List<NetworkEdge> particleState,
                          List<LineagePair> coalescibleLineagePairs) {
        // Select lineages to coalesce

        LineagePair coalescingPair = coalescibleLineagePairs.get(
                Randomizer.nextInt(coalescibleLineagePairs.size()));


        // Create coalescent node
        NetworkNode coalescentNode = new NetworkNode();
        coalescentNode.setHeight(coalescentTime)
                .addChildEdge(coalescingPair.lineage1)
                .addChildEdge(coalescingPair.lineage2);
        coalescingPair.lineage1.parentNode = coalescentNode;
        coalescingPair.lineage2.parentNode = coalescentNode;

        // Merge segment flags:
        BitSet hasSegments = new BitSet();
        hasSegments.or(coalescingPair.lineage1.hasSegments);
        hasSegments.or(coalescingPair.lineage2.hasSegments);

        // Create new lineage
        NetworkEdge lineage = new NetworkEdge(null, coalescentNode, hasSegments);
        coalescentNode.addParentEdge(lineage);

        particleState.remove(coalescingPair.lineage1);
        particleState.remove(coalescingPair.lineage2);
        particleState.add(lineage);
    }

    private void reassort(double reassortmentTime, List<NetworkEdge> particleState) {
        NetworkEdge lineage = particleState.get(Randomizer.nextInt(particleState.size()));

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
        NetworkNode node = new NetworkNode();
        node.setHeight(reassortmentTime).addChildEdge(lineage);

        // Create reassortment lineages
        NetworkEdge leftLineage = new NetworkEdge(null, node, hasSegs_left);
        NetworkEdge rightLineage = new NetworkEdge(null, node, hasSegs_right);
        node.addParentEdge(leftLineage);
        node.addParentEdge(rightLineage);

        particleState.remove(lineage);
        particleState.add(leftLineage);
        particleState.add(rightLineage);
    }


    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Not implemented.");
    }
}
