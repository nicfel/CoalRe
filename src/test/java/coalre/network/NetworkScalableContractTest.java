package coalre.network;

import org.junit.Assert;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;

/**
 * Asserts the BEAST 3 {@code Scalable} contract for {@link Network}.
 *
 * <p>The contract (beast.base.inference.Scalable) requires three mutually
 * consistent operations on a single dilation axis. For Network the axis is the
 * sum of internal-node margins, matching its interval-scaling {@code scale}.
 */
public class NetworkScalableContractTest {

    private static final double EPS = 1e-9;

    /** Two taxa with one reassortment node (#1), so both scale branches are exercised. */
    private static final String NEWICK =
            "(((A[&segments={0,1,2,3,4,5,6,7}]:1)#1[&segments={1,4}]:1,"
            + "B[&segments={0,1,2,3,4,5,6,7}]:2)[&segments={0,1,2,3,4,5,6,7}]:1,"
            + "#1[&segments={0,2,3,5,6,7}]:2)[&segments={0,1,2,3,4,5,6,7}]:0.0;";

    private Network network() {
        return new Network(NEWICK);
    }

    private static Map<String, Double> leafHeights(Network n) {
        Map<String, Double> heights = new HashMap<>();
        for (NetworkNode node : n.getNodes())
            if (node.isLeaf())
                heights.put(node.getTaxonLabel(), node.getHeight());
        return heights;
    }

    /** Invariant 1: after scale(s), getScalableValue() is exactly s x its previous value. */
    @Test
    public void scaleEquivariance() {
        for (double s : new double[]{0.25, 0.5, 1.0, 1.5, 3.0}) {
            Network n = network();
            double v0 = n.getScalableValue();
            n.scale(s);
            Assert.assertEquals("sum of margins must scale by exactly s (s=" + s + ")",
                    s * v0, n.getScalableValue(), EPS * Math.max(1.0, s * v0));
        }
    }

    /** Invariant 2: setScalableValue(V) is a fixed point of getScalableValue(). */
    @Test
    public void setIsFixedPointOfGet() {
        for (double V : new double[]{0.5, 2.0, 7.25}) {
            Network n = network();
            n.setScalableValue(V);
            Assert.assertEquals(V, n.getScalableValue(), EPS * Math.max(1.0, V));
        }
    }

    /** Invariant 3: setScalableValue(get() * s) lands in the same state as scale(s). */
    @Test
    public void setComposesWithScale() {
        double s = 1.7;

        Network viaScale = network();
        viaScale.scale(s);

        Network viaSet = network();
        viaSet.setScalableValue(viaSet.getScalableValue() * s);

        Assert.assertEquals(viaScale.toString(), viaSet.toString());
    }

    /** scale() returns the log Jacobian determinant, NOT a degrees-of-freedom count. */
    @Test
    public void scaleReturnsLogJacobian() {
        double s = 2.0;
        Network n = network();
        long dof = n.getInternalNodes().stream().filter(x -> !x.isLeaf()).count();

        double logJ = n.scale(s);

        Assert.assertEquals("log Jacobian must be dof * log(s)",
                dof * Math.log(s), logJ, EPS);
        Assert.assertNotEquals("must not return a bare dof count", (double) dof, logJ, EPS);
    }

    /** Interval scaling keeps sampled (leaf) heights fixed. */
    @Test
    public void leafHeightsPreserved() {
        Network n = network();
        Map<String, Double> before = leafHeights(n);

        n.scale(2.5);

        Map<String, Double> after = leafHeights(n);
        Assert.assertEquals(before.keySet(), after.keySet());
        for (String taxon : before.keySet())
            Assert.assertEquals("leaf " + taxon + " must not move",
                    before.get(taxon), after.get(taxon), EPS);
    }

    /** Interval scaling is valid for any positive s, so it never throws. */
    @Test
    public void neverThrowsForPositiveS() {
        for (double s : new double[]{1e-6, 0.1, 1.0, 10.0, 1e6}) {
            try {
                network().scale(s);
            } catch (RuntimeException ex) {
                Assert.fail("scale(" + s + ") must not throw, but threw " + ex);
            }
        }
    }
}
