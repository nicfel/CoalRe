package coalre.dynamics;

import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class PiecewiseConstantReassortmentRatesTest {

    private PiecewiseConstantReassortmentRates reassortmentRates;

    @Before
    public void setUp() {
        reassortmentRates = new PiecewiseConstantReassortmentRates();

        // Set up reassortment rates and rate shifts
        RealParameter reassortmentRate = new RealParameter(new Double[]{0.1, 0.2, 0.3, 0.4});
        RealParameter rateShifts = new RealParameter(new Double[]{0.0, 1.0, 2.0, 3.0});

        reassortmentRates.reassortmentRateInput.setValue(reassortmentRate, reassortmentRates);
        reassortmentRates.rateShiftsInput.setValue(rateShifts, reassortmentRates);

        reassortmentRates.initAndValidate();
    }

    @Test
    public void testGetPopSize() {
        // Test different times to ensure the correct reassortment rate is returned
        assertEquals(0.1, reassortmentRates.getPopSize(0.5), 1e-12);
        assertEquals(0.2, reassortmentRates.getPopSize(1.5), 1e-12);
        assertEquals(0.3, reassortmentRates.getPopSize(2.5), 1e-12);
        assertEquals(0.4, reassortmentRates.getPopSize(3.5), 1e-12);
    }

    @Test
    public void testGetIntegral() {
        // Test the integral over different intervals
        assertEquals(0.1 * 0.5, reassortmentRates.getIntegral(0.0, 0.5), 1e-12);
        assertEquals(0.1 * 1.0 + 0.2 * 0.5, reassortmentRates.getIntegral(0.0, 1.5), 1e-12);
        assertEquals(0.2 * 1.0 + 0.3 * 0.5, reassortmentRates.getIntegral(1.0, 2.5), 1e-12);
        assertEquals(0.1 * 1.0 + 0.2 * 1.0 + 0.3 * 1.0, reassortmentRates.getIntegral(0.0, 3.0), 1e-12);
        assertEquals(0.1 * 1.0 + 0.2 * 1.0 + 0.3 * 1.0 + 0.4 * 1.0, reassortmentRates.getIntegral(0.0, 4.0), 1e-12);

        assertEquals(0.2 * 0.5, reassortmentRates.getIntegral(1.0, 1.5), 1e-12);
        assertEquals(0.3 * 0.5, reassortmentRates.getIntegral(2.0, 2.5), 1e-12);
        assertEquals(0.4 * 0.5, reassortmentRates.getIntegral(3.0, 3.5), 1e-12);
        
        // cross individual boundaries
        assertEquals(0.1 * 0.5 + 0.2 * 0.5, reassortmentRates.getIntegral(0.5, 1.5), 1e-12);
        assertEquals(0.2 * 0.5 + 0.3 * 0.5, reassortmentRates.getIntegral(1.5, 2.5), 1e-12);
        assertEquals(0.3 * 0.5 + 0.4 * 0.5, reassortmentRates.getIntegral(2.5, 3.5), 1e-12);
        assertEquals(0.1 * 0.5 + 0.2 * 1 + 0.3 * 0.5, reassortmentRates.getIntegral(0.5, 2.5), 1e-12);
        assertEquals(0.2 * 0.5 + 0.3 * 1 + 0.4 * 0.5, reassortmentRates.getIntegral(1.5, 3.5), 1e-12);
        assertEquals(0.1 * 0.5 + 0.2 * 1 + 0.3 * 1 + 0.4 * 0.5, reassortmentRates.getIntegral(0.5, 3.5), 1e-12);
        assertEquals(0.4 * 1 , reassortmentRates.getIntegral(4.5, 5.5), 1e-12);        
        assertEquals(0.4 * 2 , reassortmentRates.getIntegral(3.5, 5.5), 1e-12);

    }
}
