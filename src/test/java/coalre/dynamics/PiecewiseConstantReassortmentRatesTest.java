package coalre.dynamics;

import org.junit.Assert;
import org.junit.Test;

import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealVectorParam;

public class PiecewiseConstantReassortmentRatesTest {

    private PiecewiseConstantReassortmentRates reassortmentRates;
    private RealVectorParam<Real> reassortmentRateInput;
    private RealVectorParam<NonNegativeReal> rateShiftsMock;

    @Test
    public void testPiecwiseConstant() {
        reassortmentRates = new PiecewiseConstantReassortmentRates();
        reassortmentRateInput = new RealVectorParam<>(new double[]{0, -1, -2, -3}, Real.INSTANCE);
        rateShiftsMock = new RealVectorParam<>(new double[]{0, 1, 2, 3}, NonNegativeReal.INSTANCE);
        reassortmentRates.initByName("reassortmentRate", reassortmentRateInput, "rateShifts", rateShiftsMock);
        
        Assert.assertEquals(Math.exp(-2) + Math.exp(-3), reassortmentRates.getIntegral(2.0, 4.0), 1e-10);        
        Assert.assertEquals(Math.exp(-1) * 0.5 + Math.exp(-2) + Math.exp(-3)*7, reassortmentRates.getIntegral(1.5, 10.0), 1e-10);
        
        
        Assert.assertEquals(Math.exp(-1) * 0.5 + Math.exp(-2) + Math.exp(-3)*7, reassortmentRates.getIntegral(1.5, 10.0), 1e-10);
        
    }
}
