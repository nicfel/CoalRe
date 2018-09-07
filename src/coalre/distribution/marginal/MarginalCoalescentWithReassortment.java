package coalre.distribution.marginal;

import beast.core.Distribution;
import beast.core.Function;
import beast.core.Input;
import beast.core.State;
import beast.util.Randomizer;

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
    double[] logParticleWeights;


    public MarginalCoalescentWithReassortment() { }

    @Override
    public void initAndValidate() {
        eventList = eventListInput.get();
        nParticles = nParticlesInput.get();

        reassortmentRate = reassortmentRateInput.get();
        populationSize = populationSizeInput.get();

        particleStates = new ParticleState[nParticles];
        particleStatesPrime = new ParticleState[nParticles];
        logParticleWeights = new double[nParticles];
        for (int p=0; p<nParticles; p++) {
            particleStates[p] = new ParticleState(eventList.getNSegments());
            particleStatesPrime[p] = new ParticleState(eventList.getNSegments());
            logParticleWeights[p] = 0.0;
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

            double maxLogWeight = Double.NEGATIVE_INFINITY;
            for (int p=0; p<nParticles; p++) {

                logParticleWeights[p] = particleStates[p].propagateParticleState(
                        prevEventTime, event.time, reassortmentRate, populationSize);

                if (logParticleWeights[p]>Double.NEGATIVE_INFINITY)
                    logParticleWeights[p] += particleStates[p].getObservedEventProbability(event, populationSize);

                if (logParticleWeights[p] > maxLogWeight)
                    maxLogWeight = logParticleWeights[p];
            }

            if (maxLogWeight == Double.NEGATIVE_INFINITY) {
                logP = Double.NEGATIVE_INFINITY;
                break;
            }

            // Compute average weight

            double sumScaledWeights = 0.0;
            for (int p=0; p<nParticles; p++) {
                logParticleWeights[p] = logParticleWeights[p] - maxLogWeight;
                sumScaledWeights += Math.exp(logParticleWeights[p]);
            }

            logP += Math.log(sumScaledWeights/nParticles) + maxLogWeight;

            // Resample particles

            double[] u = new double[nParticles];
            for (int p=0; p<nParticles; p++)
                u[p] = Randomizer.nextDouble()*sumScaledWeights;

            Arrays.sort(u);

            double accumulator = 0;
            int pPrime = 0;
            for (int p=0; p<nParticles; p++) {

                while (pPrime<nParticles && (u[pPrime]-accumulator < Math.exp(logParticleWeights[p]))) {
                    particleStates[p].copyParticleState(particleStatesPrime[pPrime]);
                    pPrime += 1;
                }

                accumulator += Math.exp(logParticleWeights[p]);
            }

            ParticleState[] tmp = particleStates;
            particleStates = particleStatesPrime;
            particleStatesPrime = tmp;

            prevEventTime = event.time;
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
