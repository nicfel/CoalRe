package coalre.dynamics;

import org.junit.Assert;
import org.junit.Test;

import beast.base.inference.parameter.RealParameter;

public class PiecewiseConstantReassortmentRateFromSkygrowthTest {

    @Test
    public void testConstantRates() {
        // Create the test instance
    	SkygrowthReassortmentRatesFromSkygrowthNe sky = new SkygrowthReassortmentRatesFromSkygrowthNe();

        // These parameters make the (log) population size plus mapping constant across intervals:
        //   logNe = 0, neToReassortment = 0 => rate = 0
        //   Because growth depends on differences between consecutive rates,
        //   growth = 0 in each interval, giving a piecewise-constant function of exp(0) = 1.
        RealParameter logNe = new RealParameter("1 -1 2");
        RealParameter neToReassortment = new RealParameter("5 -2 3");
        RealParameter rateShifts = new RealParameter("0 1 2"); // Two breakpoints => three intervals

        // Initialize the sky object
        sky.initByName(
            "logNe", logNe, 
            "neToReassortment", neToReassortment, 
            "rateShifts", rateShifts
        );

        // Check population size within each relevant interval:
        Assert.assertEquals(Math.exp(6.0), sky.getPopSize(0.0), 1e-12);
        Assert.assertEquals(Math.exp(-3.0), sky.getPopSize(1.0), 1e-12);
        Assert.assertEquals(Math.exp(1.5), sky.getPopSize(0.5), 1e-12);
        Assert.assertEquals(Math.exp(1.0), sky.getPopSize(1.5), 1e-12);
        Assert.assertEquals(Math.exp(5.0), sky.getPopSize(10.0), 1e-12);
        

        // Verify integrals of the population function:
        double integral1 = Math.exp(6.0)/-9 * Math.exp(-9.0*1.0) - Math.exp(6)/-9.0 * Math.exp(-9.0*0.0);
        Assert.assertEquals(integral1, sky.getIntegral(0.0, 1), 1e-12);

        double integral2 = Math.exp(-3)/8 * Math.exp(8.0*0.5) - Math.exp(-3)/8.0 * Math.exp(8.0*0.1);
        Assert.assertEquals(integral2, sky.getIntegral(1.1, 1.5), 1e-12);
        Assert.assertEquals(Math.exp(5.0)*7, sky.getIntegral(3.0, 10.0), 1e-12);
        
        // test inverse intensity
        for (int i = 0; i < 100; i++) {
            double t = Math.random() * 10;
            double integral = sky.getIntegral(0, t);
            double invIntensity = sky.getInverseIntensity(integral);
            Assert.assertEquals(t, invIntensity, 1e-12);
        
        }
       
        
    }
}
