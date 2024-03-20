package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.linear.*;

import java.io.PrintStream;


/**
 * @author Nicola F. Mueller
 */
@Description("Populaiton function with values at certain time points that are interpolated in between. Parameter has to be in log space")
public class Spline extends CalculationNode implements Loggable {
    final public Input<RealParameter> InfectedInput = new Input<>("logInfected",
            "Nes over time in log space", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts",
            "When to switch between elements of Ne", Input.Validate.REQUIRED);
    final public Input<RealParameter> uninfectiousRateInput = new Input<>("uninfectiousRate",
            "Rate at which individuals become uninfectious", Input.Validate.REQUIRED);
    final public Input<Integer> gridPointsInput = new Input<>("gridPoints",
            "Number of grid points to use for the spline calculation", 1000);
    final public Input<Boolean> infectedIsNeInput = new Input<>("infectedIsNe",
            "Whether the infected parameter is actually the number of infected or the logNe", false);

    RealParameter infected;
    RealParameter rateShifts;
    RealParameter uninfectiousRate;
    int gridPoints;

    double[] transmissionRate;
    double[] transmissionRateStored;

    double[] I;
    double[] I_stored;

    double[][] splineCoeffs;
    double[][] splineCoeffs_stored;

    double[] time;

    boolean ratesKnows=false;
    boolean isValid = true;
    boolean infectedIsNe = false;

    @Override
    public void initAndValidate() {
        infected = InfectedInput.get();
        rateShifts = rateShiftsInput.get();
        infected.setDimension(rateShifts.getDimension());
        uninfectiousRate = uninfectiousRateInput.get();
        gridPoints = gridPointsInput.get();
        infectedIsNe = infectedIsNeInput.get();
        recalculateRates();
    }

    // computes the Ne's at the break points from the growth rates and the transmission rates
    private void recalculateRates() {
//        notAKnotCubicSpline();
        computeAkimaSplineCoefficients();

//        clampedCubicSpline();
        // make the time grid from 0 to rateShifts.getArrayValue(rateShifts.getDimension()-1) using gridPoints
        time = new double[gridPoints+1];
        I = new double[gridPoints+1];
        transmissionRate = new double[gridPoints+1];
        double dt = rateShifts.getArrayValue(rateShifts.getDimension()-1) / (time.length-1);
        int j = 0;
        int k = j-1;
        isValid = true;
        for (int i=0; i < gridPoints; i++) {
            // update the time for this grid point
            time[i] = i*dt;
            // find the interval in which this grid point lies
            if (time[i] >= rateShifts.getArrayValue(j)) {
                j++;
                k++;
                if (k==rateShifts.getDimension()-1) {
                    k--;
                }
            }
            // get the time diff from the last point where logI was estimated
            double timeDiff = time[i]-rateShifts.getArrayValue(k);
            double timeDiff2 = timeDiff*timeDiff;
            double timeDiff3 = timeDiff2*timeDiff;
            // compute the number of infected individuals at the grid points
            I[i] = Math.exp(splineCoeffs[k][0]*timeDiff3 + splineCoeffs[k][1]*timeDiff2 + splineCoeffs[k][2]*timeDiff + splineCoeffs[k][3]);
            // compute the transmission rate at the grid points, from dI/dt and the recovery rate, the minus in front
            // of I is because the transmission rates are forward in time, but the dI/dt is backward in time
            if (infectedIsNe) {
                transmissionRate[i] = 1;
            }else {
                transmissionRate[i] = uninfectiousRate.getValue() -
                        (3 * splineCoeffs[k][0] * timeDiff2 + 2 * splineCoeffs[k][1] * timeDiff + splineCoeffs[k][2]);
            }
            if (transmissionRate[i] < 0) {
//                System.out.println(time[i] + " " + I[i] + " " + transmissionRate[i] + " " +  splineCoeffs[k][0]*timeDiff3 + " " + splineCoeffs[k][1]*timeDiff2 + " " + splineCoeffs[k][2]*timeDiff + " " + splineCoeffs[k][3]);
                isValid = false;
            }
        }
        // set the last time point to the last time point in the spline
        time[gridPoints] = Double.POSITIVE_INFINITY;
        // set the last I to the last I in the spline
        I[gridPoints] = I[gridPoints-1];
        // set the last transmission rate to the last transmission rate in the spline
        transmissionRate[gridPoints] = transmissionRate[gridPoints-1];
        ratesKnows = true;
    }

    public boolean update() {
        if (!ratesKnows) {
            recalculateRates();
        }
        return isValid;
    }

    @Override
    public boolean requiresRecalculation() {
        ratesKnows = false;
        return true;
    }

    @Override
    public void store() {
        transmissionRateStored = new double[transmissionRate.length];
        System.arraycopy(transmissionRate, 0, transmissionRateStored, 0, transmissionRate.length);
        I_stored = new double[I.length];
        System.arraycopy(I, 0, I_stored, 0, I.length);
        super.store();
    }

    @Override
    public void restore() {
        recalculateRates();
        super.restore();
    }



    /** computes the coefficients for the cubic spline interpolation
     *  de Boor, Carl. A Practical Guide to Splines. Springer-Verlag, New York: 1978
     */
    public void notAKnotCubicSpline() {
        int n = rateShifts.getDimension();

        // Calculate h values (difference between x values)
        double[] h = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            h[i] = rateShifts.getArrayValue(i + 1) - rateShifts.getArrayValue(i);
        }

        // Calculate the difference in y values
        double[] delta = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            delta[i] = (infected.getArrayValue(i + 1) - infected.getArrayValue(i)) / h[i];
        }

        // Create the tridiagonal system
        RealMatrix A = new Array2DRowRealMatrix(n, n);
        RealVector r = new ArrayRealVector(n);

        // Set up the system A*mu = r for the not-a-knot condition
        A.setEntry(0, 0, h[1]);
        A.setEntry(0, 1, -(h[0] + h[1]));
        A.setEntry(0, 2, h[0]);
        A.setEntry(n - 1, n - 3, h[n - 2]);
        A.setEntry(n - 1, n - 2, -(h[n - 2] + h[n - 3]));
        A.setEntry(n - 1, n - 1, h[n - 3]);

        for (int i = 1; i < n - 1; i++) {
            A.setEntry(i, i - 1, h[i - 1]);
            A.setEntry(i, i, 2 * (h[i - 1] + h[i]));
            A.setEntry(i, i + 1, h[i]);
            r.setEntry(i, 6 * (delta[i] - delta[i - 1]));
        }

        // Solve for mu
        DecompositionSolver solver = new LUDecomposition(A).getSolver();
        RealVector mu = solver.solve(r);

        // Store coefficients for each spline segment
        splineCoeffs = new double[n - 1][4];
        for (int i = 0; i < n - 1; i++) {
            splineCoeffs[i][0] = (mu.getEntry(i + 1) - mu.getEntry(i)) / (6 * h[i]);
            splineCoeffs[i][1] = mu.getEntry(i) / 2;
            splineCoeffs[i][2] = delta[i] - h[i] * (2 * mu.getEntry(i) + mu.getEntry(i + 1)) / 6;
            splineCoeffs[i][3] = infected.getArrayValue(i);
        }
    }

    /**
     * Clamped cubic spline interpolation with specified first derivatives at the end points.
     */
    public void clampedCubicSpline() {
        int n = rateShifts.getDimension();

        // Calculate h values (difference between x values)
        double[] h = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            h[i] = rateShifts.getArrayValue(i + 1) - rateShifts.getArrayValue(i);
        }

        // Calculate the difference in y values
        double[] delta = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            delta[i] = (infected.getArrayValue(i + 1) - infected.getArrayValue(i)) / h[i];
        }

        // Create the tridiagonal system
        RealMatrix A = new Array2DRowRealMatrix(n, n);
        RealVector r = new ArrayRealVector(n);

        // get the derivate of the first intervals
        double ddtstart = (infected.getArrayValue(1)-infected.getArrayValue(0))/(rateShifts.getArrayValue(1)-rateShifts.getArrayValue(0));
        // get the derivate of the last intervals
        double ddtend = (infected.getArrayValue(n-1)-infected.getArrayValue(n-2))/(rateShifts.getArrayValue(n-1)-rateShifts.getArrayValue(n-2));

        // Clamped condition at the start
        A.setEntry(0, 0, 2 * h[0]);
        A.setEntry(0, 1, h[0]);
        r.setEntry(0, ((delta[0] - ddtstart) * 3));

        // Clamped condition at the end
        A.setEntry(n - 1, n - 2, h[n - 2]);
        A.setEntry(n - 1, n - 1, 2 * h[n - 2]);
        r.setEntry(n - 1, ((ddtend - delta[n - 2]) * 3));

        // Set up the middle equations of the tridiagonal system
        for (int i = 1; i < n - 1; i++) {
            A.setEntry(i, i - 1, h[i - 1]);
            A.setEntry(i, i, 2 * (h[i - 1] + h[i]));
            A.setEntry(i, i + 1, h[i]);
            r.setEntry(i, 3 * (delta[i] - delta[i - 1]));
        }

        // Solve for mu
        DecompositionSolver solver = new LUDecomposition(A).getSolver();
        RealVector mu = solver.solve(r);

        // Store coefficients for each spline segment
        splineCoeffs = new double[n - 1][4];
        for (int i = 0; i < n - 1; i++) {
            splineCoeffs[i][0] = (mu.getEntry(i + 1) - mu.getEntry(i)) / (6 * h[i]);
            splineCoeffs[i][1] = mu.getEntry(i) / 2;
            splineCoeffs[i][2] = delta[i] - h[i] * (2 * mu.getEntry(i) + mu.getEntry(i + 1)) / 6;
            splineCoeffs[i][3] = infected.getArrayValue(i);
        }
    }


    private void computeAkimaSplineCoefficients() {
        int n = rateShifts.getDimension();
        double[] slopes = new double[n-1];

        // Compute central differences for slopes
        for (int i = 0; i < (n-1); i++) {
            slopes[i] = (infected.getArrayValue(i + 1) - infected.getArrayValue(i )) / (rateShifts.getArrayValue(i + 1) - rateShifts.getArrayValue(i));
        }

        // Calculate the Akima weights and the weighted slopes
        double[] weightedSlopes = new double[n];
        for (int i = 2; i < (n - 2); i++) {
            double weight1 = Math.abs(slopes[i + 1] - slopes[i]);
            double weight2 = Math.abs(slopes[i-1] - slopes[i - 2]);
            if (weight1 + weight2 == 0)
                weightedSlopes[i] = (slopes[i + 1] + slopes[i]) / 2;
            else
                weightedSlopes[i] = (weight1 * slopes[i-1] + weight2 * slopes[i]) / (weight1 + weight2);
        }

        weightedSlopes[0] = slopes[0];
        weightedSlopes[1] = (slopes[0] + slopes[1])/2;

        // Extrapolate the end slopes
        weightedSlopes[n-1] = slopes[n-2];
        weightedSlopes[n - 2] = (slopes[n-3]-slopes[n-2])/2;

        // Calculate coefficients for Akima spline segments
        splineCoeffs = new double[n - 1][4];
        for (int i = 0; i < n - 1; i++) {
            splineCoeffs[i] = computeAkimaCoefficients(rateShifts.getArrayValue(i), rateShifts.getArrayValue(i + 1),
                    infected.getArrayValue(i), infected.getArrayValue(i + 1),
                    weightedSlopes[i], weightedSlopes[i + 1]);
        }
    }

    /**
     * Method to compute coefficients for an Akima spline segment.
     */
    private double[] computeAkimaCoefficients(double x0, double x1, double y0, double y1, double d0, double d1) {
        double dx = x1 - x0;
        double[] coeffs = new double[4];
        coeffs[3] = y0;
        coeffs[2] = d0;
        coeffs[1] = (3 * (y1 - y0) / dx - 2 * d0 - d1) / dx;
        coeffs[0] = (d0 + d1 - 2 * (y1 - y0) / dx) / (dx * dx);
        return coeffs;
    }

    @Override
    public void init(PrintStream out) {
        for (int i = 0; i < splineCoeffs.length; i++){
            for (int j = 0; j < splineCoeffs[i].length; j++){
                out.print("splineCoeffs_" + i + "_" + j + "\t");
            }
        }
    }

    @Override
    public void log(long l, PrintStream printStream) {
        for (int i = 0; i < splineCoeffs.length; i++){
            for (int j = 0; j < splineCoeffs[i].length; j++){
                printStream.print(splineCoeffs[i][j] + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream printStream){

    }



    }