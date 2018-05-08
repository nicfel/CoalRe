package coalre.distribution;

import java.util.ArrayList;
import java.util.List;


import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import coalre.network.NetworkIntervalType;
import coalre.network.NetworkNode;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a beast.tree using under the framework of Mueller (2017).")
public class CoalescentWithReassortment extends NetworkDistribution {
	
	public Input<RealParameter> rRateInput = new Input<>("reassortmentRate", "reassortment rate", Input.Validate.REQUIRED);
	public Input<RealParameter> coalescentRateInput = new Input<>("coalescentRate", "coalescent rate", Input.Validate.REQUIRED);
	

    // Set up for lineage state probabilities
    private ArrayList<NetworkNode> activeLineageNodes;
	private ArrayList<Boolean> activeLineageIsLeft;
    private ArrayList<Boolean[]> activeSegments;
    
    private double coalescentRate;
    private double reassortmentRate;
           
    @Override
    public void initAndValidate(){    	
    }
    
//    int count = 0;
    public double calculateLogP() {
//    	System.out.println("calc");
    	logP = 0;
    	// newly calculate tree intervals
    	List<NetworkEvent> networkEventList = networkIntervalsInput.get().getNetworkEventList();
    	double nextEventTime = 0.0;

    	boolean first = true;
    	
    	activeLineageNodes = new ArrayList<>();
    	activeLineageIsLeft = new ArrayList<>();
    	activeSegments = new ArrayList<>();
    	
    	coalescentRate = coalescentRateInput.get().getValue(0);
    	reassortmentRate = rRateInput.get().getValue(0);

    	for (NetworkEvent event : networkEventList) {
        	nextEventTime = event.time;
        	
        	if (nextEventTime > 0) {
        		logP += intervalContribution(nextEventTime);
        	}

        	switch (event.type) {
				case COALESCENCE:
					logP += coalesce(event);
					break;

				case SAMPLE:
					sample(event, first);
					break;

				case REASSORTMENT:
					logP += reassortmentstart(event);
					reassortmentend(event);
					break;
			}
			first = false;

       		if (logP==Double.NEGATIVE_INFINITY)
       			break;
        }
        
		return logP;
    }
    
	private double reassortmentstart(NetworkEvent event) {
		List<NetworkLineage> reassLineages = event.lineagesRemoved;

    	if (reassLineages.size() != 1) {
			System.err.println("Unsupported number of incoming lineages at a reassortment event");
			System.exit(0);
		}
    	
    	final int daughterIndex = activeLineages.indexOf(reassLineages.get(0).getNr());
    	
		if (daughterIndex == -1 ) {
			System.out.println(reassLineages.get(0).getNr() + " " + activeLineages);
			System.out.println("daughter lineage at reassortment event not found");
			return Double.NaN;
		}


    	Boolean[] segsBefore = reassLineages.get(0).getHasSegments();
    	
    	int c_before = 0;
    	for (int i = 0; i < segsBefore.length; i++){
    		if (segsBefore[i]==true) c_before++;
    	}
    	
    	
        activeLineages.add(reassLineages.get(0).getParent().getNr());
        
        activeSegments.add(reassLineages.get(0).getParent().getHasSegments());
        
        activeLineages.remove(daughterIndex);
        activeSegments.remove(daughterIndex);


		return Math.log(reassortmentRate*2*Math.pow(0.5, c_before));
		
	}


	private void reassortmentend(NetworkEvent event) {
		List<NetworkLineage> reassLineages = event.lineagesAdded;
		
    	if (reassLineages.size() != 1) {
			System.err.println("Unsupported number of incoming lineages at a second part of a reassortment event");
			System.exit(0);
    	}
		
		activeLineages.add(reassLineages.get(0).getNr());
		activeSegments.add(reassLineages.get(0).getHasSegments());
	}


	private void sample(NetworkEvent event, boolean first) {
		List<NetworkLineage> incommingLineages = event.lineagesAdded;
		
		for (NetworkNode l : incommingLineages) {
			activeLineages.add(l.getNr());
			activeSegments.add(l.getHasSegments());
		}
	}

	private double coalesce(NetworkEvent event) {
		List<NetworkNode> coalLines = networkIntervalsInput.get().getLineagesRemoved(networkInterval);
		
    	if (coalLines.size() > 2) {
			System.err.println("Unsupported coalescent at non-binary node");
			System.exit(0);
		}
    	if (coalLines.size() < 2) {
    		System.out.println();
    		System.out.println("WARNING: Less than two lineages found at coalescent event!");
    		System.out.println();
    		return Double.NaN;
		}
		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println(coalLines.get(0).getNr() + " " + coalLines.get(1).getNr() + " " + activeLineages);
			System.out.println("daughter lineages at coalescent event not found");
			return Double.NaN;
		}
		
        activeLineages.add(coalLines.get(0).getParent().getNr());    
        
        activeSegments.add(coalLines.get(0).getParent().getHasSegments());
        
		activeLineages.remove(Math.max(daughterIndex1, daughterIndex2));
		activeLineages.remove(Math.min(daughterIndex1, daughterIndex2));
		
		activeSegments.remove(Math.max(daughterIndex1, daughterIndex2));
		activeSegments.remove(Math.min(daughterIndex1, daughterIndex2));

		
		return Math.log(coalescentRate);
	}

	private double intervalContribution(double nextEventTime) {
		// for each active lineage, calculate the reassortment interval contribution
		double reassCont = 0.0;
		// TODO make this thing more efficient
		for (int i = 0; i < activeSegments.size(); i++){
	    	int nr_lins = 0;
	    	for (int j = 0; j < activeSegments.get(i).length; j++)
	    		if (activeSegments.get(i)[j]==true) nr_lins++;
	    	
	    	// Precompute powers
	    	reassCont -= (1-2*Math.pow(0.5, nr_lins));

		}
		reassCont *= reassortmentRate;
		
		
		int overlap = 0;
		
		// TODO keep directly track of which network nodes overlap
		for (int i = 0; i < activeSegments.size(); i++){
			for (int j = i+1; j < activeSegments.size(); j++){
				for (int k = 0; k < activeSegments.get(i).length; k++){
					if (activeSegments.get(i)[k] && activeSegments.get(i)[k]){
						overlap++;
						break; // break out of the loop if there is at least one segment overlapping						
					}// if
				}// k					
			}// j		
		}// i
		
		double coalCont = -overlap*coalescentRate;
		
		return nextEventTime*(reassCont+coalCont);
	}   

    
}
