package coalre.simulator;

import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

public class SISwithReassortment extends Network implements Loggable {
    public Input<Double> transmissionRateInput = new Input<>("transmissionRate", "transmission rate", Input.Validate.REQUIRED);

    public Input<Double> recoveryRateInput  = new Input<>("recoveryRate", "recovery rate", Input.Validate.REQUIRED);

    public Input<Double> reassortmenProbabilityInput = new Input<>("reassortmenProbability", "reassortment probability", Input.Validate.REQUIRED);

    public Input<Double> samplingProbabilityInput = new Input<>("samplingProbability", "sampling probability", Input.Validate.REQUIRED);

    public Input<Integer> populationSizeInput = new Input<>("populationSize", "population size", Input.Validate.REQUIRED);

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
    int nSegments;

    int states;

    List<SIREvents> lineages;
    List<SIREvents> sirEvents;

    int S,I,R;

    @Override
    public void initAndValidate(){

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
        List<Individual> sampledIndividuals = new ArrayList<>();

        // set the starting conditions for simulations
        do {
            sirEvents = new ArrayList<>();

            double S = populationSizeInput.get() - 1;
            double I = 1;
            double R = 0;
            int individualNr = 0;
            // start the simulations while keeping track of the different individuals
            List<Individual> activeIndividuals = new ArrayList<>();
            sampledIndividuals = new ArrayList<>();

            Individual root = new Individual(individualNr);
            individualNr++;
            activeIndividuals.add(root);
            double time = 0.0;
            do {
                double transmissionRate = I*(S + I - 1)/(S+I) * transmissionRateInput.get() * Math.max(0,(S-R)/S);
                double recoveryRate = I * recoveryRateInput.get();
                

                double totalRate = transmissionRate + recoveryRate;
                double nextEventTime = Math.log(1.0 / Randomizer.nextDouble()) / totalRate;

                // pick which event happens next based on the rates, by drawing at random
                double random = Randomizer.nextDouble() * totalRate;

                if (random < transmissionRate) {
                    // transmission event
                    // pick a random individual to transmit to
                    Individual parent = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()));
                    // set time for parent
                    parent.setTime(time + nextEventTime);
                    // create two new individuals
                    Individual child1 = new Individual(individualNr);
                    individualNr++;
                    Individual child2 = new Individual(individualNr);
                    individualNr++;
                    // start the first new lineage
                    parent.addChild(child1);
                    child1.addParent(parent);
                    // add the child to the second parent
                    parent.addChild(child2);
                    child2.addParent(parent);

                    // add the child to the list of active individuals
                    activeIndividuals.add(child1);

                    // pick if the individual to be infected
                    double probCoInf = S / (double) (S + I - 1);
//                    System.out.println(probCoInf);
                    if (Randomizer.nextDouble() > probCoInf) {
                        // co-infection event
                        // pick another active individual that is not the parent and not child1
                        Individual parent2 = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()-1));
                        while (parent2 == parent) {
                            parent2 = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()-1));
                        }
                        parent2.setTime(time + nextEventTime);

                        // add the child to the second parent
                        child2.setTime(time + nextEventTime);

                        // create a third child
                        Individual child3 = new Individual(individualNr);
                        individualNr++;
                        child3.addParent(child2);
                        child2.addChild(child3);

                        child3.addParent(parent2);
                        parent2.addChild(child3);

                        activeIndividuals.add(child3);
                        activeIndividuals.remove(parent2);
                        sirEvents.add(new SIREvents(0, time + nextEventTime));
                    } else {
                        // update the population counts
                        S--;
                        I++;
                        activeIndividuals.add(child2);
                        sirEvents.add(new SIREvents(1, time + nextEventTime));
                    }
                    // remove the parent
                    activeIndividuals.remove(parent);
                } else if (random < (transmissionRate + recoveryRate)) {
                    // pick a random individual to recover
                    Individual individual = activeIndividuals.get(Randomizer.nextInt(activeIndividuals.size()));
                    individual.setTime(time + nextEventTime);
                    // choose if that individual will be sampled
                    if (Randomizer.nextDouble() < samplingProbabilityInput.get()) {
                        // sample the individual
                        sampledIndividuals.add(individual);
                    }
                    // remove the individual from the list of active individuals
                    activeIndividuals.remove(individual);
                    // update the population counts
                    I--;
                    R++;
                    S++;
                    sirEvents.add(new SIREvents(2, time + nextEventTime));
                } else {
                    System.out.println(transmissionRate + recoveryRate);
                    System.out.println(transmissionRate + " " + recoveryRate);
                    System.out.println(I);
                    System.out.println((S + I) );
                    System.exit(0);
                    // waning immunity event
                    // update the population counts
                    R--;
                    S++;
                    sirEvents.add(new SIREvents(3, time + nextEventTime));
                }
                time += nextEventTime;
//                System.out.println("S: " + S + " I: " + I + " R: " + R + " " + activeIndividuals.size() + " " + sampledIndividuals.size());
            } while (time < simulationTimeInput.get() && I > 0);
            System.out.println("start building network from " + sampledIndividuals.size() + " sampled individuals" + " simulation time was " + time);
            System.out.println("total number of recovered was " + R );

            // build the network from the sampled individuals
        }while (sampledIndividuals.size()<minSamplesInput.get());
        buildNetwork(sampledIndividuals);
//        System.out.println(lineages);
    }

    protected void buildNetwork(List<Individual> sampledIndividuals) {
        // get the time of the most recent sample
        double mrsi_time = 0.0;
        for (Individual individual:sampledIndividuals){
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

        List<Individual> activeIndividuals = new ArrayList<>();
        List<NetworkEdge> activeEdges = new ArrayList<>();
        System.out.print("start");
        while (activeIndividuals.size()>1 || sampledIndividuals.size()>0){
//            System.out.println("get next sampling time");
            // get the next sampling event time, while keeping track of the index of the individual
            double nextSamplingTime = Double.NEGATIVE_INFINITY;
            int nextSamplingIndex = -1;
            for (int i=0;i<sampledIndividuals.size();i++){
                if(sampledIndividuals.get(i).getTime()>nextSamplingTime){
                    nextSamplingTime = sampledIndividuals.get(i).getTime();
                    nextSamplingIndex = i;
                }
            }
//            System.out.println("get next active time");
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
//                    else if ((activeIndividuals.get(i).getParents().get(0).getTime() == nextActiveTime)){
//                        // if the other individual has two parents, then do nothing
//                        if (activeIndividuals.get(nextActiveIndex).getParents().size()==2){
//                            System.out.println("do nothing");}
//                        else{
//                            System.out.println(activeIndividuals.get(nextActiveIndex).getParents().size());
//                            System.out.println(activeIndividuals.get(i).getParents().size());
//                            System.out.println("----------------");
//                        }
//
//                    }
                }
            }
//            System.out.println(" done " + nextSamplingTime +" " + nextActiveTime + "");


            // depending on which one is next
            if (nextSamplingTime>nextActiveTime){
//                System.out.println("sampling event");
                NetworkNode sampledNode = new NetworkNode();
                sampledNode.setHeight(mrsi_time-nextSamplingTime);
                sampledNode.setTaxonLabel("sample_no"+sampledIndividualNr);
                sampledNode.setTaxonIndex(sampledIndividualNr);
                sampledIndividualNr++;

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

                lineages.add(new SIREvents(0, nextSamplingTime));
            }else {
                // check if the individual has two parents
                if (activeIndividuals.get(nextActiveIndex).getParents().size()==2){
//                    System.out.println("co-infection event");
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

//                        System.out.println("observed reassortment event");
                        // make a reassortment node
                        NetworkNode reassortmentNode = new NetworkNode();
                        reassortmentNode.setHeight(mrsi_time-nextActiveTime);
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
                        lineages.add(new SIREvents(1, nextActiveTime, prob));
                    }else{
                        // the event is not observed, just replace the individual with its parent that has cardinality>0
                        if (segmentsLeft.cardinality()>0) {
                            activeIndividuals.set(nextActiveIndex, activeIndividuals.get(nextActiveIndex).getParents().get(0));
                        }else{
                            activeIndividuals.set(nextActiveIndex, activeIndividuals.get(nextActiveIndex).getParents().get(1));
                        }
                    }
                }else{
                    // get the other child for this coalescent event
                    Individual otherChild = activeIndividuals.get(nextActiveIndex).getParents().get(0).getChildren().get(0);
                    if (otherChild.equals(activeIndividuals.get(nextActiveIndex))) {
                        otherChild = activeIndividuals.get(nextActiveIndex).getParents().get(0).getChildren().get(1);
                    }
                    // check if the other child is in the activeIndividuals list
                    if (activeIndividuals.contains(otherChild)) {
                        totalObservedCoal++;
                        // make a coalescent node
                        NetworkNode coalNode = new NetworkNode();
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
                        lineages.add(new SIREvents(2, nextActiveTime));
                    } else {
                        // replace the activeIndividuals with its parent
                        activeIndividuals.set(nextActiveIndex, activeIndividuals.get(nextActiveIndex).getParents().get(0));
                    }
                }
            }
            // check if all the individuals in activeIndividuals have different numbers
            for (int i=0; i<activeIndividuals.size();i++){
                for (int j=i+1; j<activeIndividuals.size();j++){
                    if (activeIndividuals.get(i).equals(activeIndividuals.get(j))){
                        System.out.println(activeIndividuals);
                        System.out.println(activeIndividuals.get(i));
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


    protected class Individual{
        List<Individual> parent;
        List<Individual> child;

        double time;
        int number;


        int S, I;

        public Individual(int number){
            parent = new ArrayList<>();
            child = new ArrayList<>();
            this.number = number;
        }

        public List<Individual> getParents(){
            return parent;
        }

        public List<Individual> getChildren(){
            return child;
        }

        public double getTime(){
            return time;
        }

        public void setTime(double time){
            this.time = time;
        }

        public void addParent(Individual parent){
            this.parent.add(parent);
        }

        public void addChild(Individual child){
            this.child.add(child);
        }

        // build a method to check if two individual objects ar the same
        public boolean equals(Individual other){
        	if (this.number == other.number){
        		return true;
        	}else{
        		return false;
        	}
        }

        // define what happens in toString
        public String toString(){
        	return "" + number ;
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

    protected class SIREvents{
        int eventType;
        double time;

        double prob = 0.0;

        public SIREvents(int eventType, double time){
            this.eventType = eventType;
            this.time = time;
        }

        public SIREvents(int eventType, double time, double prob){
            this.eventType = eventType;
            this.time = time;
            this.prob = prob;
        }

        @Override
        public String toString() {
            return eventType + ":" + time + ":" + prob;
        }
    }

    @Override
    public void init(PrintStream out) {
        out.println("SIR\tlineages\t");
    }

    @Override
    public void close(PrintStream out) {

    }

    @Override
    public void log(long sample, PrintStream out) {
        out.println(sirEvents + "\t" + lineages+"\t");
    }
}

