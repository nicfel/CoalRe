package coalre.dynamics;
import coalre.dynamics.PiecewiseConstantReassortmentRateScalers;
import beast.base.inference.parameter.RealParameter;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class PiecewiseConstantReassortmentRateScalersTest {

    private PiecewiseConstantReassortmentRateScalers scalers;
    private Spline spline;
    private RealParameter infectedToRho;
    private RealParameter rateShifts;

    @Before
    public void setUp() {
        spline = new Spline();
        
        
        
        
        infectedToRho = new RealParameter(new Double[]{1.0, 2.0});
        rateShifts = new RealParameter(new Double[]{0.0, 1.0});
        
        scalers = new PiecewiseConstantReassortmentRateScalers();
        scalers.splineInput.setValue(spline, scalers);
        scalers.InfectedToRhoInput.setValue(infectedToRho, scalers);
        scalers.rateShiftsInput.setValue(rateShifts, scalers);
        
        scalers.initAndValidate();
    }

    @Test
    public void testGetPopSize() {
        assertEquals(10.0, scalers.getPopSize(0.0), 1e-6);
        assertEquals(20.0, scalers.getPopSize(1.0), 1e-6);
        assertEquals(60.0, scalers.getPopSize(2.0), 1e-6);
    }

    @Test
    public void testGetIntegral() {
        double integral = scalers.getIntegral(0.0, 2.0);
        assertTrue(integral > 0);
        // You would add specific checks based on the expected integral calculation.
    }

    @Test
    public void testGetIntensity() {
        double intensity = scalers.getIntensity(2.0);
        assertTrue(intensity > 0);
        // You would add specific checks based on the expected intensity calculation.
    }
}
