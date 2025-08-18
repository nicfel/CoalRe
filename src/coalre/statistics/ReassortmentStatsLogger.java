package coalre.statistics;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Tree;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

public class ReassortmentStatsLogger extends BEASTObject implements Loggable {


    public Input<Network> networkInput = new Input<>("network",
            "Network for which to log statistics.",
            Input.Validate.REQUIRED);
    

    Network network;
    int segCount;
    boolean logObservable = false;


    public ReassortmentStatsLogger() { }

    @Override
    public void initAndValidate() {
        network = networkInput.get();
        segCount = network.getSegmentCount();

    }

    @Override
    public void init(PrintStream out) {
        String prefix = network.getID() == null ? "networkStat." : network.getID() + ".";

//		if (logObservable){
	        for (int i = 0; i < segCount; i++){
//	        	for (int j = i+1; j < segCount; j++){
                out.print(prefix + "segmentTreeHeight." + i + "\t");        		
//	                out.print(prefix + "splitReassortment." + i + "_" + j + "\t");        		
//	        	}
	        }    
//		}else{
//	        for (int i = 0; i < segCount; i++){
//	        	for (int j = i+1; j < segCount; j++){
//	                out.print(prefix + "jointObsReassortment." + i + "_" + j + "\t");        		
//	                out.print(prefix + "splitObsReassortment." + i + "_" + j + "\t");        		
//	        	}
//	        }    
//
//		}
		
    }

    @Override
    public void log(long sample, PrintStream out) {
		for (int i = 0; i < segCount; i++) {
			out.print(NetworkStatsLogger.getHeightSegmentsRoot(network.getRootEdge(), i) + "\t");
		}
    	
  

	}

    @Override
    public void close(PrintStream out) {

    }
}
