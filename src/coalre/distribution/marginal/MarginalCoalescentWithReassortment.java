package coalre.distribution.marginal;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.*;

public class MarginalCoalescentWithReassortment extends Distribution {

    public Input<RealParameter> reassortmentRateInput = new Input<>(
	        "reassortmentRate",
            "reassortment rate (per lineage per unit time)",
            Input.Validate.REQUIRED);

	public Input<PopulationFunction> populationFunctionInput = new Input<>(
	        "populationModel",
            "Population model.",
            Input.Validate.REQUIRED);

	public Input<ObservedEventList> eventListInput = new Input<>(
	        "observedEventList",
            "List of observed segment tree events.",
            Input.Validate.REQUIRED);

    public Input<Integer> nParticlesInput = new Input<>(
            "nParticles",
            "Number of particles to use in SMC calculation.",
            Input.Validate.REQUIRED);

    RealParameter reassortmentRate;
    PopulationFunction populationFunction;

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
        populationFunction = populationFunctionInput.get();

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
        double logP = 0.0;

        double currentTime = startTime;
        while (true) {
            int k = particleState.size();

            double currentTransformedTime = populationFunction.getIntensity(currentTime);
            double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
            double timeToNextCoal = populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;

            double timeToNextReass = k>=1 ? Randomizer.nextExponential(k*reassortmentRate.getValue()) : Double.POSITIVE_INFINITY;

            double nextEventTime = Math.min(timeToNextCoal, timeToNextReass);

            if (nextEventTime>nextObservationTime)
                break;


        }

        return logP;
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
