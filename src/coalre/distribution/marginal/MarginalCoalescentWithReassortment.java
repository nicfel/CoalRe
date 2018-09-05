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
    ParticleState[] particleStates, particleStatesPrime;
    double[] particleWeights;


    public MarginalCoalescentWithReassortment() { }

    @Override
    public void initAndValidate() {
        eventList = eventListInput.get();
        nParticles = nParticlesInput.get();

        reassortmentRate = reassortmentRateInput.get();
        populationSize = populationSizeInput.get();

        particleStates = new ParticleState[nParticles];
        for (int p=0; p<nParticles; p++) {
            particleStates[p] = new ParticleState(eventList.getNSegments());
            particleStatesPrime[p] = new ParticleState(eventList.getNSegments());
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

                particleWeights[p] = particleStates[p].propagateParticleState(
                        prevEventTime, event.time, reassortmentRate, populationSize);

                if (particleWeights[p]>Double.NEGATIVE_INFINITY)
                    particleWeights[p] += particleStates[p].getObservedEventProbability(event,
                            reassortmentRate, populationSize);

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
                    particleStates[p].copyParticleState(particleStatesPrime[pPrime]);
                    pPrime += 1;
                }

                accumulator += Math.exp(particleWeights[p]);
            }

            particleStates = particleStatesPrime;

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

    @Override
    public boolean isStochastic() {
        return true;
    }
}
