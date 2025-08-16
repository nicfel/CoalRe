package coalre.simulator;

import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

public class SuperspreadingStructuredSIRwithReassortment extends Network implements Loggable {

    public Input<RealParameter> transmissionRateInput = new Input<>("transmissionRate", "transmission rate", Input.Validate.REQUIRED);

    public Input<RealParameter> recoveryRateInput  = new Input<>("recoveryRate", "recovery rate", Input.Validate.REQUIRED);

    public Input<RealParameter> waningImmunityRateInput  = new Input<>("waningImmunityRate", "waning immunity rate", Input.Validate.REQUIRED);

    public Input<RealParameter> reassortmenProbabilityInput = new Input<>("reassortmenProbability", "reassortment probability", Input.Validate.REQUIRED);

    public Input<RealParameter> samplingProbabilityInput = new Input<>("samplingProbability", "sampling probability", Input.Validate.REQUIRED);

    public Input<IntegerParameter> populationSizeInput = new Input<>("populationSize", "population size", Input.Validate.REQUIRED);

    public Input<Integer> nSegmentsInput = new Input<>("nSegments","Number of segments. Used if no segment trees are supplied.");

    public Input<Double> simulationTimeInput = new Input<>("simulationTime", "simulation time", Input.Validate.REQUIRED);

    public Input<String> fileNameInput = new Input<>("fileName",
            "Name of file to write simulated network to.");

    public Input<Boolean> enableSegTreeUpdateInput = new Input<>("enableSegmentTreeUpdate",
            "If false, segment tree objects won't be updated to agree with simulated " +
                    "network. (Default true.)", true);

    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "One or more segment trees to initialize.", new ArrayList<>());

    public Input<Integer> minSamplesInput = new Input<>("minSamples",
            "Minimum number of samples in the tree, otherwise, simulation is repeated.", 0);
    public Input<Double> kInput = new Input<>("k",
            "k-value for the negative binomial distribution", Input.Validate.REQUIRED);

    public Input<RealParameter> migrationRatesInput = new Input<>("migrationRates",
            "Migration rates between segments", Input.Validate.REQUIRED);

    int nSegments;
    int states;

    List<StructuredSIREvents> lineages;

    List<StructuredSIREvents> structuredSirEvents;

    int I[];
    int S[];
    int R[];

    public void init(PrintStream out) {
        out.println("SIR\tlineages\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.println(structuredSirEvents + "\t" + lineages+"\t");
    }

    @Override
    public void close(PrintStream out) {

    }



    @Override
    public void initAndValidate(){
        states = transmissionRateInput.get().getDimension();

        List<String> IDs = new ArrayList<>();
        if (segmentTreesInput.get().isEmpty()) {
            nSegments = (int) nSegmentsInput.get();
        }else {
            nSegments = segmentTreesInput.get().size();
            segmentNames = new String[nSegments];
            // initialize names of segments
            baseName = "";
            for (int segIdx=0; segIdx<nSegments; segIdx++) {
                IDs.add(segmentTreesInput.get().get(segIdx).getID());
                segmentNames[segIdx] = segmentTreesInput.get().get(segIdx).getID().replace(baseName, "");
            }
//            segmentTreesInput.get().clear();
        }

        simulateNetwork();
        // Create the taxon sets from the sampled individuals
        List<Taxon> taxonList = new ArrayList<>();
        for (NetworkNode leaf : getLeafNodes()) {
            Taxon taxon = new Taxon(leaf.getTaxonLabel());
            taxonList.add(taxon);
        }
        // make a TaxonSet from the taxa
        TaxonSet taxonset = new TaxonSet();
        taxonset.initByName("taxon", taxonList);

        // create the trait set using the sampling times
        TraitSet traitSet = new TraitSet();
        // make a string of values with leaf.getTaxonLabel()+"=" +leaf.getHeight() and individual entries seprated by ,
        String values = "";
        for (NetworkNode leaf : getLeafNodes()) {
            values += leaf.getTaxonLabel()+"=" +leaf.getHeight() + ",";
        }
        traitSet.initByName("value", values.substring(0, values.length()-1), "taxa", taxonset, "traitname", "date-backward");

        // Update segment trees:
        if (enableSegTreeUpdateInput.get()) {
        	System.out.println(this);
            for (int segIdx = 0; segIdx < nSegments; segIdx++) {
            	
                // get the segment tree, set the taxonsets and update the segment tree based on the network
                Tree segmentTree = new Tree();
                segmentTree.initByName("taxonset", taxonset);
                updateSegmentTree(segmentTree, segIdx);
                segmentTree.setID(IDs.get(segIdx));
//                segmentTreesInput.get().add( segmentTree);
                segmentTree.setEverythingDirty(false);
                segmentTreesInput.get().set(segIdx, segmentTree);
            }
        }
        // Write simulated network to file if requested
        super.initAndValidate();

        for (int segIdx = 0; segIdx < nSegments; segIdx++) {
            System.out.println(segmentTreesInput.get().get(segIdx).getID());
            System.out.println(segmentTreesInput.get().get(segIdx).getRoot().toNewick());
            System.out.println(segmentTreesInput.get().get(segIdx).getNodeCount());

        }

        // Write simulated network to file if requested
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {

                ps.println(toString().replace("0.0)", "0.0000000001)"));
                for (int segIdx = 0; segIdx < nSegments; segIdx++) {
                    ps.println(segmentTreesInput.get().get(segIdx).getRoot().toNewick() +";");
                }
            } catch (FileNotFoundException ex) {
                throw new RuntimeException("Error writing to output file '"
                        + fileNameInput.get() + "'.");
            }
        }
        super.initAndValidate();
    }

    protected void simulateNetwork() {
        List<StructuredIndividual> sampledIndividuals = new ArrayList<>();

        // set the starting conditions for simulations
        do {
            structuredSirEvents = new ArrayList<>();
            // sample initial location
            int initLoc = Randomizer.nextInt(states);
            S = new int[states];
            I = new int[states];
            R = new int[states];
            for (int i = 0; i < states; i++) {
                if (i == initLoc) {
                    S[i] = populationSizeInput.get().getValue(i) - 1;
                    I[i] = 1;
                    R[i] = 0;
                } else {
                    S[i] = populationSizeInput.get().getValue(i);
                    I[i] = 0;
                    R[i] = 0;
                }
            }
            int individualNr = 0;
            // start the simulations while keeping track of the different individuals
            List<StructuredIndividual> activeIndividuals = new ArrayList<>();
            sampledIndividuals = new ArrayList<>();

            StructuredIndividual root = new StructuredIndividual(individualNr, initLoc);
            individualNr++;
            activeIndividuals.add(root);
            double time = 0.0;
            System.out.println("start");
            do {
                double transmissionRate = 0;
                double recoveryRate = 0;
                double waningRate = 0;
                for (int i = 0; i < states; i++) {
                    transmissionRate += I[i] * recoveryRateInput.get().getValue(i); // This should be recovery rate, as we model the Reff through the avg number of secondary infections below
                    recoveryRate += I[i] * recoveryRateInput.get().getValue(i);
                    waningRate += R[i] * waningImmunityRateInput.get().getValue(i);
                }
                double migrationRate = 0;
                int c = 0;
                for (int i = 0; i < states; i++) {
                    for (int j = 0; j < states; j++) {
                        if (i!=j){
                            migrationRate += I[i] * migrationRatesInput.get().getValue(c);
                            c++;
                        }
                    }
                }

                double totalRate = transmissionRate + recoveryRate + waningRate + migrationRate;
                double nextEventTime = Math.log(1.0 / Randomizer.nextDouble()) / totalRate;

                // pick which event happens next based on the rates, by drawing at random
                double random = Randomizer.nextDouble() * totalRate;

                if (random < transmissionRate) {
//                    System.out.println("t");
                    ReturnVal rv = transmit(activeIndividuals, transmissionRate, time, nextEventTime, individualNr);
                    time = rv.time;
                    individualNr = rv.individualNr;
                } else if (random < (transmissionRate + recoveryRate)) {
//                    System.out.println("r");
                    individualNr = recover(activeIndividuals, recoveryRate, time, nextEventTime, individualNr, sampledIndividuals);
                }else if (random < (transmissionRate + recoveryRate + migrationRate)){
//                    System.out.println("m");
                    individualNr = migrate(activeIndividuals, migrationRate, time, nextEventTime, individualNr);
                } else {
                    System.out.println(transmissionRate + recoveryRate);
                    System.out.println(transmissionRate + " " + recoveryRate + " " + waningRate);
                    System.out.println(I);
//                    System.out.println((S + I) );
                    System.exit(0);
                    // waning immunity event
                    // update the population counts
//                    R--;
//                    S++;
//                    structuredSirEvents.add(new SIREvents(3, time + nextEventTime));
                }
                time += nextEventTime;
                // log every 1000th iteration
                if (individualNr % 1000 == 0)
                    System.out.println(time + " S: " + getSumI(S) + " I: " + getSumI(I) + " R: " + getSumI(R) + " samples: " + sampledIndividuals.size() + " active: " + activeIndividuals.size());
            } while (time < simulationTimeInput.get() && getSumI(I) > 0);
            System.out.println("start building network from " + sampledIndividuals.size() + " sampled individuals" + " simulation time was " + time);
            // build the network from the sampled individuals
        }while (sampledIndividuals.size()<minSamplesInput.get());
        
        buildNetwork(sampledIndividuals);
//        System.out.println(lineages);
    }

    private int getSumI(int[] I){
        int sum = 0;
        for (int i = 0; i<states; i++){
            sum+=I[i];
        }
        return sum;
    }

    private ReturnVal transmit(List<StructuredIndividual> activeIndividuals, double transmissionRate,
                          double time, double nextEventTime, int individualNr){
        // pick the location of the transmission event
        double randomT = Randomizer.nextDouble();
        double cummulative = 0;
        for (int i = 0; i < states; i++){
            cummulative += I[i] * recoveryRateInput.get().getValue(i)/transmissionRate;
            if (randomT<=cummulative){
                // pick the number of offsprings from a negative binomial distribution with R and k
                // from a gamma and a poisson distribution
                double secondary_infections = transmissionRateInput.get().getArrayValue(i)/recoveryRateInput.get().getArrayValue(i);
                double gamma = Randomizer.nextGamma(kInput.get(), kInput.get() / (double) secondary_infections);
                int nOffspring = (int) Randomizer.nextPoisson(gamma);
                int isR = 0;
                // for each offspring, randomly sample if it hits an R
                for (int j = 0; j < nOffspring; j++) {
                    if (Randomizer.nextDouble() < (R[i]+1) / (double) (S[i] + I[i] + R[i])) {
                        // sample the individual
                        isR++;
                    }
                }
                nOffspring -= isR;
                System.out.println();

                // don't do anything if there are no offspring
                if (nOffspring > 0){
                    // transmission event
                    // pick a random individual to transmit to
//                    System.out.println("a " + activeIndividuals + " " + i + " " + Arrays.toString(I) + " " + transmissionRate);
                    StructuredIndividual parent = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()));
                    while (parent.type!=i){
                        parent = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()));
                    }

                    // build the superspreading event as a series of regular infection events
                    for (int j = 0; j < nOffspring; j++) {
                        // set time for parent
                        parent.setTime(time + nextEventTime);
                        // create two new individuals
                        StructuredIndividual child1 = new StructuredIndividual(individualNr, i);
                        individualNr++;
                        StructuredIndividual child2 = new StructuredIndividual(individualNr, i);
                        individualNr++;
                        // start the first new lineage
                        parent.addChild(child1);
                        child1.addParent(parent);
                        // add the child to the second parent
                        parent.addChild(child2);
                        child2.addParent(parent);
                        System.out.println(parent);

                        // add the child to the list of active individuals
                        activeIndividuals.add(child1);

                        // pick if the individual to be infected, assuming an individual is co-infected once at most
                        double probCoInf = S[i] / (double) (S[i] + I[i]-1-j);
                        if (Randomizer.nextDouble() > probCoInf) {
                            // offset for the number of active individuals to draw from
                            int totOptions = j==0? 1 : 2 + j;

                            // co-infection event
                            // pick another active individual that is not the parent, while making sure it is not
                            // the same as any of the previous individuals in this event

                            StructuredIndividual parent2 = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()-totOptions));
                            while (parent2 == parent || parent2.type!=i) {
                                parent2 = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()-totOptions));
                            }
//                            System.out.println(i + " " + parent2);

                            parent2.setTime(time + nextEventTime);

                            // add the child to the second parent
                            child2.setTime(time + nextEventTime);

                            // create a third child
                            StructuredIndividual child3 = new StructuredIndividual(individualNr, i);

                            individualNr++;
                            child3.addParent(child2);
                            child2.addChild(child3);

                            child3.addParent(parent2);
                            parent2.addChild(child3);

                            activeIndividuals.add(child3);
                            activeIndividuals.remove(parent2);
                            
                    		System.out.println("co-inf" + parent2 + " " + child3);
                            
                            structuredSirEvents.add(new StructuredSIREvents(0, time + nextEventTime, i));
                        } else {
//                            System.out.println("...");

                            // update the population counts
                            S[i]--;
                            I[i]++;
                            activeIndividuals.add(child2);
                            structuredSirEvents.add(new StructuredSIREvents(1, time + nextEventTime, i));
//                            System.out.println("...");

                        }
                        // check that the last two active individuals are both in type i
						if (activeIndividuals.get(activeIndividuals.size() - 1).type != i
								|| activeIndividuals.get(activeIndividuals.size() - 2).type != i) {
							System.out.println("error");
							System.exit(0);
						}
						System.out.println(child1 + " " + child2 + " " + parent);
                        	
                        // remove the parent
                        activeIndividuals.remove(parent);
                        parent = child1;
                        time += 0.0000000001;
                    }
                }
				if (Math.abs(time - 64.22849979972362+nextEventTime) < 0.00001) {
					System.out.println(activeIndividuals);
					// also print the parents of all
					System.out.print("[");
					for (StructuredIndividual individual : activeIndividuals) {
						System.out.print(individual.getParents().get(0) + ", ");
					}
//					System.exit(0);
				}
                
                
                return new ReturnVal(time, individualNr);
            }
        }

        return new ReturnVal(time, individualNr);
    }

    private int migrate(List<StructuredIndividual> activeIndividuals,
                         double migrationRate, double time, double nextEventTime, int individualNr){
        // pick the route of migration
        double randomMig = Randomizer.nextDouble();
        double cummulative = 0;
        int c = 0;
        for (int i = 0; i<states; i++){
            for (int j = 0; j < states; j++){
                if (i!=j){
                    cummulative+= (I[i] * migrationRatesInput.get().getValue(c))/migrationRate;
                    if (randomMig<cummulative){
                        StructuredIndividual individual = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()));
                        while (individual.type!=i){
                            individual = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()));
                        }
                        StructuredIndividual child = new StructuredIndividual(individualNr, j);
                        individualNr++;
                        individual.addChild(child);
                        child.addParent(individual);
                        individual.setTime(time + nextEventTime);
                        activeIndividuals.add(child);
                        activeIndividuals.remove(individual);
                        structuredSirEvents.add(new StructuredSIREvents(3, time + nextEventTime, new int[]{i,j}));

                        I[i]--;
                        I[j]++;

                        return individualNr;
                    }
                    c++;
                }
            }
        }
        return individualNr;

    }

    private int recover(List<StructuredIndividual> activeIndividuals, double recoveryRate,
                         double time, double nextEventTime, int individualNr,
                         List<StructuredIndividual> sampledIndividuals){
        // pick the state to recover from
        double randomRec = Randomizer.nextDouble();
        double cummulative = 0;
        for (int i = 0; i<states; i++){
            cummulative+= (I[i] * recoveryRateInput.get().getValue(i))/recoveryRate;
            if (randomRec<=cummulative){
                StructuredIndividual individual = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()));
//                System.out.println(i + " " + activeIndividuals + " " + Arrays.toString(I));
                while (individual.type!=i){
                    individual = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()));
                }
                individual.setTime(time + nextEventTime);
                // choose if that individual will be sampled
                if (Randomizer.nextDouble() < samplingProbabilityInput.get().getArrayValue(i)) {
                    // sample the individual
                    sampledIndividuals.add(individual);
                }
                activeIndividuals.remove(individual);
                I[i]--;
                R[i]++;
                structuredSirEvents.add(new StructuredSIREvents(2, time + nextEventTime, i));
                return individualNr;
            }
        }
        return individualNr;
    }
    
    protected void buildNetwork(List<StructuredIndividual> sampledIndividuals) {
        // get the time of the most recent sample
        double mrsi_time = 0.0;
        for (StructuredIndividual individual : sampledIndividuals){
            if(individual.getTime()>mrsi_time){
                mrsi_time = individual.getTime();
            }
        }

        System.out.println("mrsi time: " + mrsi_time);
        System.out.println("=================");
        System.out.println("=================");
        System.out.println("=================");
        System.out.println("=================");
        System.out.println("=================");
        System.out.println("=================");
        
        lineages = new ArrayList<>();
        

        int sampledIndividualNr = 0;
        int totalObservedCoal=0;

        List<StructuredIndividual> activeIndividuals = new ArrayList<>();
        List<NetworkEdge> activeEdges = new ArrayList<>();
        System.out.println("start");
        while (activeIndividuals.size()>1 || sampledIndividuals.size()>0){

            // get the next sampling event time, while keeping track of the index of the individual
            double nextSamplingTime = Double.NEGATIVE_INFINITY;
            int nextSamplingIndex = -1;
            for (int i=0;i<sampledIndividuals.size();i++){
                if(sampledIndividuals.get(i).getTime()>nextSamplingTime){
                    nextSamplingTime = sampledIndividuals.get(i).getTime();
                    nextSamplingIndex = i;
                }
            }

            // get the activeIndividuals time, while keeping track of the index of the individual
            double nextActiveTime = Double.NEGATIVE_INFINITY;
            int nextActiveIndex = -1;
            for (int i=0;i<activeIndividuals.size();i++){
                if (activeIndividuals.get(i).getParents().size()>0) {
                    if (activeIndividuals.get(i).getParents().size() == 2){
                        double nextParentTime = Math.max(activeIndividuals.get(i).getParents().get(0).getTime(),
                                activeIndividuals.get(i).getParents().get(1).getTime());
                        if (nextParentTime>=nextActiveTime) {
                            nextActiveTime = nextParentTime;
                            nextActiveIndex = i;
                        }
                    } else if (activeIndividuals.get(i).getParents().get(0).getTime() > nextActiveTime) {
                        nextActiveTime = activeIndividuals.get(i).getParents().get(0).getTime();
                        nextActiveIndex = i;
                    }
                }
            }
            
            // depending on which one is next
            if (nextSamplingTime>nextActiveTime){
                NetworkNode sampledNode = new NetworkNode();
                sampledNode.setHeight(mrsi_time-nextSamplingTime);
                sampledNode.setTaxonLabel("sample_no"+sampledIndividualNr);
                sampledNode.setTaxonIndex(sampledIndividualNr);
                sampledNode.setMetaData(",type="+sampledIndividuals.get(nextSamplingIndex).type);
                sampledIndividualNr++;

                lineages.add(new StructuredSIREvents(0, nextSamplingTime, sampledIndividuals.get(nextSamplingIndex).type));
                // set segments for the sampled node, where all bits are true
                BitSet segments = new BitSet(nSegments);
                // set all bits to true
                segments.set(0,nSegments);
                // make a new edge
                NetworkEdge edge = new NetworkEdge(null, sampledNode, segments);
                sampledNode.addParentEdge(edge);
                activeEdges.add(edge);
                activeIndividuals.add(sampledIndividuals.get(nextSamplingIndex));

                // remove the individual from the list of sampled individuals
                sampledIndividuals.remove(nextSamplingIndex);

            }else {
//                 check if the individual has two parents
                if (activeIndividuals.get(nextActiveIndex).getParents().size()==2){
                    System.out.println("co-infection event");
                    // get the corresponding edge
                    NetworkEdge edge = activeEdges.get(nextActiveIndex);
                    // pick randomly which segments go left and which go right
                    BitSet segmentsLeft = new BitSet(nSegments);
                    BitSet segmentsRight = new BitSet(nSegments);
                    for (int i=0;i<nSegments;i++){
                        if (edge.hasSegments.get(i)){
                            if (Randomizer.nextBoolean()){
                                segmentsLeft.set(i);
                            }else{
                                segmentsRight.set(i);
                            }
                        }
                    }
                    // check if the two segments have cardinality > 0
                    if (segmentsLeft.cardinality()>0 && segmentsRight.cardinality()>0){
                        // compute the average observation probability
                        double prob = 0.0;
                        for (NetworkEdge e : activeEdges){
                            prob += 1-2*Math.pow(0.5, e.hasSegments.cardinality());
                        }
                        prob = prob/activeEdges.size();

                        System.out.println("observed reassortment event");
                        // make a reassortment node
                        NetworkNode reassortmentNode = new NetworkNode();
                        reassortmentNode.setMetaData(",type="+activeIndividuals.get(nextActiveIndex).type);
                        reassortmentNode.setHeight(mrsi_time-nextActiveTime);
                        lineages.add(new StructuredSIREvents(1, nextActiveTime, prob, activeIndividuals.get(nextActiveIndex).type));

                        // set the child edge
                        edge.parentNode = reassortmentNode;
                        reassortmentNode.addChildEdge(edge);
                        // make two new edges
                        NetworkEdge edgeLeft = new NetworkEdge(null, reassortmentNode, segmentsLeft);
                        NetworkEdge edgeRight = new NetworkEdge(null, reassortmentNode, segmentsRight);
                        // add the edges to the activeEdges list
                        reassortmentNode.addParentEdge(edgeLeft);
                        reassortmentNode.addParentEdge(edgeRight);
                        activeEdges.add(edgeLeft);
                        activeEdges.add(edgeRight);

                        edgeLeft.childNode=reassortmentNode;
                        edgeRight.childNode=reassortmentNode;

                        // add the parent individuals to the active individuals
                        activeIndividuals.add(activeIndividuals.get(nextActiveIndex).getParents().get(0));
                        activeIndividuals.add(activeIndividuals.get(nextActiveIndex).getParents().get(1));

                        // remove the individual from the activeIndividuals list
                        activeIndividuals.remove(nextActiveIndex);
                        activeEdges.remove(nextActiveIndex);
                    }else{
                        // the event is not observed, just replace the individual with its parent that has cardinality>0
                        if (segmentsLeft.cardinality()>0) {
                            activeIndividuals.set(nextActiveIndex, activeIndividuals.get(nextActiveIndex).getParents().get(0));
                        }else{
                            activeIndividuals.set(nextActiveIndex, activeIndividuals.get(nextActiveIndex).getParents().get(1));
                        }
                    }
                }else if (activeIndividuals.get(nextActiveIndex).getParents().get(0).getChildren().size()==2){
                    // get the other child for this coalescent event
                    StructuredIndividual otherChild = activeIndividuals.get(nextActiveIndex).getParents().get(0).getChildren().get(0);
                    if (otherChild.equals(activeIndividuals.get(nextActiveIndex))) {
                        otherChild = activeIndividuals.get(nextActiveIndex).getParents().get(0).getChildren().get(1);
                    }
                    // check if the other child is in the activeIndividuals list
                    if (activeIndividuals.contains(otherChild)) {
                    	System.out.println(nextActiveTime);
                        System.out.println("coalescent event "  + activeIndividuals.size() + " " + activeIndividuals);
						for (int l = 0; l < activeIndividuals.size(); l++) {
							System.out.print(activeIndividuals.get(l) + " ");
						}
						System.out.println();
						for (int l = 0; l < activeIndividuals.size(); l++) {
							System.out.print(activeIndividuals.get(l).getParents().get(0) + " ");
						}
						System.out.println();
                        totalObservedCoal++;
                        // make a coalescent node
                        NetworkNode coalNode = new NetworkNode();
                        coalNode.setMetaData(",type="+activeIndividuals.get(nextActiveIndex).type);
                        lineages.add(new StructuredSIREvents(2, nextActiveTime, activeIndividuals.get(nextActiveIndex).type));


                        coalNode.setHeight(mrsi_time - nextActiveTime);
                        // get the two child edges
                        NetworkEdge childEdge1 = activeEdges.get(nextActiveIndex);
                        NetworkEdge childEdge2 = activeEdges.get(activeIndividuals.indexOf(otherChild));
                        // add the parent node to the edges
                        childEdge1.parentNode = coalNode;
                        childEdge2.parentNode = coalNode;
                        // add the edges to the coalNode
                        coalNode.addChildEdge(childEdge1);
                        coalNode.addChildEdge(childEdge2);
                        // make a new parent edge with the segments being the union of the two child edge segments
                        BitSet segments = (BitSet) childEdge1.hasSegments.clone();
                        segments.or(childEdge2.hasSegments);
                        NetworkEdge parentEdge = new NetworkEdge(null, coalNode, segments);
                        // add the parent edge to the activeEdges list
                        activeEdges.add(parentEdge);
                        // remove the two child edges from the activeEdges list
                        activeEdges.remove(childEdge1);
                        activeEdges.remove(childEdge2);
                        activeIndividuals.add(activeIndividuals.get(nextActiveIndex).getParents().get(0));
                        // remove the two children from the activeIndividuals list
                        activeIndividuals.remove(nextActiveIndex);
                        activeIndividuals.remove(otherChild);

                    } else {
                    	// check if this is a coalescent or a migration event
                    	if (activeIndividuals.get(nextActiveIndex).type != activeIndividuals.get(nextActiveIndex).getParents().get(0).type) { 
                    		System.out.println("migration event " + activeIndividuals.size() + " " + activeIndividuals);
                    		System.out.println(activeIndividuals.get(nextActiveIndex) + " " + activeIndividuals.get(nextActiveIndex).getParents());
                    		System.out.print(activeIndividuals.get(nextActiveIndex).getParents().get(0).getChildren());
                    		System.out.println(nextActiveTime);
                    		System.exit(0);
                    	}
                        // replace the activeIndividuals with its parent
                        activeIndividuals.set(nextActiveIndex, activeIndividuals.get(nextActiveIndex).getParents().get(0));
//                		System.out.println(".......... " + activeIndividuals.size() + " " + activeIndividuals);

                    }
                }else {
                	// migration event
                    System.out.println("migration event " + activeIndividuals.size() + " " + activeIndividuals);
                    // it is a migration event, add that to the lineages
            		lineages.add(new StructuredSIREvents(3, nextActiveTime, new int[]{activeIndividuals.get(nextActiveIndex).type, activeIndividuals.get(nextActiveIndex).getParents().get(0).type}));
                    activeIndividuals.set(nextActiveIndex, activeIndividuals.get(nextActiveIndex).getParents().get(0));

                }
            }
            // check if all the individuals in activeIndividuals have different numbers
            for (int i=0; i<activeIndividuals.size();i++){
                for (int j=i+1; j<activeIndividuals.size();j++){
                    if (activeIndividuals.get(i).equals(activeIndividuals.get(j))){
                        System.out.println(activeIndividuals);
                        System.out.println(activeIndividuals.get(i));
                        System.out.print("===========dsa");
                        System.exit(0);
                    }
                }
            }
        }
        System.out.print("\n\n\n" + totalObservedCoal +" " +sampledIndividualNr + "\n\n\n");
        this.setRootEdge(activeEdges.get(0));
        System.out.println(this.getLeafNodes().size() +" " + this.getInternalNodes().size());
        System.out.println(this);
        
    }

    private class ReturnVal{
        double time;
        int individualNr;

        public ReturnVal(double time, int individualNr){
            this.time = time;
            this.individualNr = individualNr;
        }


    }


    protected class StructuredIndividual{
        List<StructuredIndividual> parent;
        List<StructuredIndividual> child;

        double time;
        int number;

        int type;


        int S, I;

        public StructuredIndividual(int number, int type){
            parent = new ArrayList<>();
            child = new ArrayList<>();
            this.number = number;
            this.type = type;
        }

        public List<StructuredIndividual> getParents(){
            return parent;
        }

        public List<StructuredIndividual> getChildren(){
            return child;
        }

        public double getTime(){
            return time;
        }

        public void setTime(double time){
            this.time = time;
        }

        public void addParent(StructuredIndividual parent){
            this.parent.add(parent);
        }

        public void addChild(StructuredIndividual child){
            this.child.add(child);
        }

        // build a method to check if two individual objects ar the same
        public boolean equals(StructuredIndividual other){
            if (this.number == other.number){
                return true;
            }else{
                return false;
            }
        }

        // define what happens in toString
        public String toString(){
            return "" + number + ":" + type;
        }

        public void setS(int S){
            this.S = S;
        }
        public int getS(){
            return S;
        }

        public void setI(int I){
            this.I = I;
        }

        public int getI(){
            return I;
        }

    }

    protected class StructuredSIREvents{
        int eventType;
        double time;
        double prob = 0.0;

        int[] fromto;

        int type;

        public StructuredSIREvents(int eventType, double time, int type){
            this.eventType = eventType;
            this.time = time;
            this.type = type;
        }

        public StructuredSIREvents(int eventType, double time, double prob, int type){
            this.eventType = eventType;
            this.time = time;
            this.prob = prob;
            this.type = type;
        }

        public StructuredSIREvents(int eventType, double time, int[] fromto){
            this.eventType = eventType;
            this.time = time;
            this.fromto = fromto;
        }


        @Override
        public String toString() {
            if (fromto!=null){
                return eventType + ":" + time + ":" + prob + ":" + fromto[0] +"_" + fromto[1];
            }else{
                return eventType + ":" + time + ":" + prob + ":" + type;
            }
        }


    }


}

