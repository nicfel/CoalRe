package coalre.util;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.HeapSort;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

public class InitFromTree extends BEASTObject implements StateNodeInitialiser {

    public Input<Network> networkInput = new Input<>("network",
            "Network used to initialize segment trees.", Input.Validate.REQUIRED);

    public Input<Tree> treeInput = new Input<>("tree",
            "Segment tree to initialize.", Input.Validate.REQUIRED);
    
    Network network;
    int nsegments;
    
	@Override
	public void initAndValidate() {
		network = networkInput.get();
		nsegments = network.getSegmentCount();
		
	}

	@Override
	public void initStateNodes() {
		
        TaxonSet taxonSet = treeInput.get().getTaxonset();
		
		calculateIntervals();
		List<NetworkEdge> extantLineages = new ArrayList<NetworkEdge>();
		List<Integer> nodeNr = new ArrayList<Integer>();
		
		double currTime = 0.0;
		for (int i = 0; i < intervals.length; i++) {
			if (lineagesRemoved[i]==null) {
				if (lineagesAdded[i]!=null)					
					for (int j = 0; j < lineagesAdded[i].size(); j++)
						sample(lineagesAdded[i].get(j), extantLineages, nodeNr, taxonSet);				
			}else {
				coalesce(lineagesAdded[i].get(0), lineagesRemoved[i], extantLineages, nodeNr);				

			}
		}
		
        network.setRootEdge(extantLineages.get(0));	
	}
	
	
    private void coalesce(Node linAdded, List<Node> linsRemoved, List<NetworkEdge> extantLineages, List<Integer> nodeNr) {
        // Sample the pair of lineages that are coalescing:
        // Create coalescent node
    	
        NetworkEdge lineage1 = extantLineages.get(nodeNr.indexOf(linsRemoved.get(0).getNr()));
        NetworkEdge lineage2 = extantLineages.get(nodeNr.indexOf(linsRemoved.get(1).getNr()));
    	
        NetworkNode coalescentNode = new NetworkNode();
        coalescentNode.setHeight(linAdded.getHeight())
                .addChildEdge(lineage1)
                .addChildEdge(lineage2);        
        
        lineage1.parentNode = coalescentNode;
        lineage2.parentNode = coalescentNode;

        // Merge segment flags:
        BitSet hasSegments = new BitSet();
        hasSegments.or(lineage1.hasSegments);
        hasSegments.or(lineage2.hasSegments);

        // Create new lineage
        NetworkEdge lineage = new NetworkEdge(null, coalescentNode, hasSegments);
        coalescentNode.addParentEdge(lineage);

        extantLineages.remove(lineage1);
        extantLineages.remove(lineage2);
        extantLineages.add(lineage);
        
        nodeNr.remove(nodeNr.indexOf(linsRemoved.get(0).getNr()));
        nodeNr.remove(nodeNr.indexOf(linsRemoved.get(1).getNr()));
        
        nodeNr.add(linAdded.getNr());

        
    }

	

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(network);		
	}
	
    private void sample(Node tn, List<NetworkEdge> extantLineages, List<Integer> nodeNr, TaxonSet taxonSet) {
        // sample the network node
        NetworkNode n = new NetworkNode();
        
        int index = taxonSet.getTaxonIndex(tn.getID());

        
        
        n.setHeight(tn.getHeight());
        n.setTaxonLabel(tn.getID());
        n.setTaxonIndex(index);
        
        nodeNr.add(tn.getNr());

        		
        // Create corresponding lineage
        BitSet hasSegs = new BitSet();
        
        hasSegs.set(0, nsegments);
        
        NetworkEdge lineage = new NetworkEdge(null, n, hasSegs);
        extantLineages.add(lineage);
        n.addParentEdge(lineage);
    }

    /**
     * Recalculates all the intervals for the given beast.tree.
     */
    @SuppressWarnings("unchecked")
    protected void calculateIntervals() {
        Tree tree = treeInput.get();

        final int nodeCount = tree.getNodeCount();

        times = new double[nodeCount];
        int[] childCounts = new int[nodeCount];

        collectTimes(tree, times, childCounts);

        indices = new int[nodeCount];

        HeapSort.sort(times, indices);

        intervals = new double[nodeCount];
        lineageCounts = new int[nodeCount];
        lineagesAdded = new List[nodeCount];
        lineagesRemoved = new List[nodeCount];

        // start is the time of the first tip
        double start = times[indices[0]];
        int numLines = 0;
        int nodeNo = 0;
        intervalCount = 0;
        while (nodeNo < nodeCount) {

            int lineagesRemoved = 0;
            int lineagesAdded = 0;

            double finish = times[indices[nodeNo]];
            double next;

            do {
                final int childIndex = indices[nodeNo];
                final int childCount = childCounts[childIndex];
                // don't use nodeNo from here on in do loop
                nodeNo += 1;
                if (childCount == 0) {
                    addLineage(intervalCount, tree.getNode(childIndex));
                    lineagesAdded += 1;
                } else {
                    lineagesRemoved += (childCount - 1);

                    // record removed lineages
                    final Node parent = tree.getNode(childIndex);
                    //assert childCounts[indices[nodeNo]] == beast.tree.getChildCount(parent);
                    //for (int j = 0; j < lineagesRemoved + 1; j++) {
                    for (int j = 0; j < childCount; j++) {
                        Node child = j == 0 ? parent.getLeft() : parent.getRight();
                        removeLineage(intervalCount, child);
                    }

                    // record added lineages
                    addLineage(intervalCount, parent);
                    // no mix of removed lineages when 0 th
                    if (multifurcationLimit == 0.0) {
                        break;
                    }
                }

                if (nodeNo < nodeCount) {
                    next = times[indices[nodeNo]];
                } else break;
            } while (Math.abs(next - finish) <= multifurcationLimit);

            if (lineagesAdded > 0) {

                if (intervalCount > 0 || ((finish - start) > multifurcationLimit)) {
                    intervals[intervalCount] = finish - start;
                    lineageCounts[intervalCount] = numLines;
                    intervalCount += 1;
                }

                start = finish;
            }

            // add sample event
            numLines += lineagesAdded;

            if (lineagesRemoved > 0) {

                intervals[intervalCount] = finish - start;
                lineageCounts[intervalCount] = numLines;
                intervalCount += 1;
                start = finish;
            }
            // coalescent event
            numLines -= lineagesRemoved;
        }
    }
    
    protected static void collectTimes(Tree tree, double[] times, int[] childCounts) {
        Node[] nodes = tree.getNodesAsArray();
        for (int i = 0; i < nodes.length; i++) {
            Node node = nodes[i];
            times[i] = node.getHeight();
            childCounts[i] = node.isLeaf() ? 0 : 2;
        }
    }
    
    protected void addLineage(int interval, Node node) {
        if (lineagesAdded[interval] == null) lineagesAdded[interval] = new ArrayList<>();
        lineagesAdded[interval].add(node);
    }

    protected void removeLineage(int interval, Node node) {
        if (lineagesRemoved[interval] == null) lineagesRemoved[interval] = new ArrayList<>();
        lineagesRemoved[interval].add(node);
    }
    
    
    /**
     * The lineages in each interval (stored by node ref).
     */
    double[] intervals;
    protected List<Node>[] lineagesAdded;
    protected List<Node>[] lineagesRemoved;
    int[] indices;
    double[] times;
    int intervalCount;
    int[] lineageCounts;    
    double multifurcationLimit=-1.0;





	
}
