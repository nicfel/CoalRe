package coalre.dynamics;

import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class PiecewiseConstantReassortmentRatesTest {

    private PiecewiseConstantReassortmentRates reassortmentRates;
    private RealParameter reassortmentRateInput;
    private RealParameter rateShiftsMock;

    @Test
    public void setUp() {
        reassortmentRates = new PiecewiseConstantReassortmentRates();
        reassortmentRateInput = new RealParameter("0 -1 -2 -3");
        rateShiftsMock = new RealParameter("0 1 2 3");
        reassortmentRates.initByName("reassortmentRate", reassortmentRateInput, "rateShifts", rateShiftsMock);
        
        System.out.println(reassortmentRates.getIntegral(0.0, 2.0));
        
        
        
        assertEquals(Math.exp(-2) + Math.exp(-3), reassortmentRates.getIntegral(2.0, 4.0));        
        assertEquals(Math.exp(-1) * 0.5 + Math.exp(-2) + Math.exp(-3)*7, reassortmentRates.getIntegral(1.5, 10.0));
        
        
        assertEquals(Math.exp(-1) * 0.5 + Math.exp(-2) + Math.exp(-3)*7, reassortmentRates.getIntegral(1.5, 10.0));

        System.out.println(reassortmentRates.getPopSize(1.0));
        System.out.println(reassortmentRates.getPopSize(2.0));
        System.out.println(reassortmentRates.getPopSize(3.0));
        System.out.println(reassortmentRates.getPopSize(5.0));
        System.out.println(reassortmentRates.getPopSize(100.0));

        System.out.println(reassortmentRates.getIntegral(1.5, 10000.0));

        
    }
}
