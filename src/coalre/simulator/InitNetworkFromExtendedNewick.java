package coalre.simulator;

import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.evolution.alignment.TaxonSet;
import coalre.network.Network;
import coalre.network.NetworkNode;

import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Initialises a Network state node from an extended Newick string supplied as
 * CDATA content of the &lt;init&gt; element.  CDATA avoids XML-escaping issues
 * with the &amp;, { and } characters present in extended Newick metadata.
 *
 * Extends Network so that the inherited fromXML() is available for state
 * restoration.  The target network is populated in initAndValidate() (not only
 * in initStateNodes()) so that auto-generated SegmentTreeInitializer objects,
 * which call network.getSegmentCount() during their own initAndValidate(), see
 * a fully populated network.
 *
 * Usage in XML:
 * <pre>
 * &lt;init spec="coalre.simulator.InitNetworkFromExtendedNewick"
 *       network="@myNetwork"&gt;
 *     &lt;![CDATA[ (...extended newick...); ]]&gt;
 * &lt;/init&gt;
 * </pre>
 */
public class InitNetworkFromExtendedNewick extends Network implements StateNodeInitialiser {

    public Input<Network> networkInput = new Input<>("network",
            "Network state node to initialise from the extended Newick string.",
            Input.Validate.REQUIRED);

    // BEAST2 passes CDATA element text content as an input named "value"
    public Input<String> valueInput = new Input<>("value",
            "Extended Newick string supplied as CDATA content of the <init> element.",
            Input.Validate.REQUIRED);

    public Input<TaxonSet> taxonSetInput = new Input<>("taxonSet",
            "Taxon set used to assign leaf taxon indices on the initialised network.",
            Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        // Populate the target network here (in addition to initStateNodes) so
        // that any SegmentTreeInitializer whose initAndValidate() runs in the
        // same initBEASTObjects pass sees a non-empty network.
        populateTargetNetwork();
    }

    @Override
    public void initStateNodes() {
        populateTargetNetwork();
    }

    private void populateTargetNetwork() {
    	
        String newick = valueInput.get();
        if (newick == null || newick.isBlank())
            throw new IllegalArgumentException(
                    "InitNetworkFromExtendedNewick: no extended Newick string found. " +
                    "Supply it as CDATA content of the <init> element.");

        newick = cleanNewick(newick.strip());
        Network parsed = new Network(newick);
        Network target = networkInput.get();
        target.setRootEdge(parsed.getRootEdge());
        target.segmentCount = null;

        // Set taxon indices on leaf nodes so that updateSegmentTree() can build
        // clade BitSets correctly (getTaxonIndex() returns -1 without this step).
        TaxonSet taxonSet = taxonSetInput.get();
        for (NetworkNode leaf : target.getLeafNodes())
            leaf.setTaxonIndex(taxonSet.getTaxonIndex(leaf.getTaxonLabel()));
    }

    /**
     * Strips metadata attributes that the Network parser does not need
     * (realCoal, segsCarried, type) and removes whitespace inside segment
     * sets so the string matches the format expected by the parser,
     * e.g. [&segments={0,1}] instead of [&segments={0, 1},realCoal=false,...].
     */
    private static String cleanNewick(String newick) {
        // Remove simulation-only metadata attributes
        newick = newick.replaceAll(",\\s*realCoal=[^,\\]]+", "");
        newick = newick.replaceAll(",\\s*segsCarried=[^,\\]]+", "");
        newick = newick.replaceAll(",\\s*type=[^,\\]]+", "");

        // Remove whitespace inside segment sets: {0, 1} -> {0,1}
        Matcher m = Pattern.compile("\\{[^}]+\\}").matcher(newick);
        StringBuffer sb = new StringBuffer();
        while (m.find()) {
            m.appendReplacement(sb, m.group().replaceAll("\\s", ""));
        }
        m.appendTail(sb);
        return sb.toString();
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(networkInput.get());
    }
}
