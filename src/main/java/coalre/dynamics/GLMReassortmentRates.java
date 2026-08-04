package coalre.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.spec.domain.NonNegativeInt;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealVector;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Nicola F. Mueller
 */
@Description("Population function with defines Ne's at points in time and interpolated between them. Parameter has to be in log space. The Ne's are used to compute the transmission rates, which are used for the coalescent process. ")
public class GLMReassortmentRates extends PopulationFunction.Abstract implements Loggable {

	final public Input<List<RealVector<Real>>> predictorInput = new Input<>("predictor", "predictor parameters that are used to calculate the Ne", new ArrayList<>());
	// written to (setDimension) in initAndValidate: needs the concrete param type
	final public Input<RealVectorParam<Real>> NeToReassortmentInput = new Input<>("neToReassortment",
			"the value that maps the number of infected or the Ne to the reassortment rate ");
	final public Input<RealVector<NonNegativeReal>> rateShiftsInput = new Input<>("rateShifts",
			"When to switch between elements of Ne", Input.Validate.REQUIRED);
	// written to (setValue) in initAndValidate: needs the concrete param type
	final public Input<IntScalarParam<NonNegativeInt>> predictorIsActiveInput = new Input<>("predictorIsActive",
			"index of the active predictor, or the number of predictors for none");
	final public Input<Integer> independentAfterInput = new Input<>("independentAfter",
			"ignore differences after that index");
	// written to (setDimension) in initAndValidate: needs the concrete param type
	final public Input<RealVectorParam<Real>> effectSizeInput = new Input<>("effectSize",
			"the effect size of the predictors on the reassortment rates", Input.Validate.REQUIRED);

	List<RealVector<Real>> predictors;
	RealVectorParam<Real> NeToReassortment;
	RealVector<NonNegativeReal> rateShifts;
	RealVectorParam<Real> effectSize;
	IntScalarParam<NonNegativeInt> predictorIsActive;

	boolean NesKnown = false;
	double[] growth;
	double[] growth_stored;

	double[] rates;
	double[] stored_rates;
	
	double[] invRatio;
	double[] invRatio_stored;
	
	@Override
	public void initAndValidate() {
		
		NeToReassortment = NeToReassortmentInput.get();
		rateShifts = rateShiftsInput.get();
		predictors = predictorInput.get();
		predictorIsActive = predictorIsActiveInput.get();
		// Valid range is [0, predictors.size()], where predictors.size() is the "no predictor
		// active" state (see the guards in recalculateNe/getStandardPredictor).
		//
		// BEAST 2 pinned this with IntegerParameter.setUpper(predictors.size()), which also
		// bounded MCMC proposals. BEAST 3 derives bounds from the Domain (plus any attached
		// distribution) rather than a settable per-instance field, so there is no equivalent
		// call. The check below catches a bad starting value, but it CANNOT constrain what
		// operators propose during the run -- that must now come from the XML, by giving
		// predictorIsActive a domain or prior that caps it at predictors.size().
		if (predictorIsActive.get() < 0 || predictorIsActive.get() > predictors.size())
			throw new IllegalArgumentException("predictorIsActive must be in [0, "
					+ predictors.size() + "] but was " + predictorIsActive.get());
		predictorIsActive.set(0);
		effectSize = effectSizeInput.get();
		effectSize.setDimension(predictors.size());

		NeToReassortment.setDimension(rateShifts.size());
		growth = new double[rateShifts.size()];

		rates = new double[rateShifts.size()];
		stored_rates = new double[rateShifts.size()];

		invRatio = new double[rateShifts.size()];
		invRatio_stored = new double[rateShifts.size()];
		
		recalculateNe();
	}

	@Override
	public List<String> getParameterIds() {
		return null;
	}

	@Override
	public double getPopSize(double t) {
		int i = getIntervalNr(t);
		double timediff = t;
		timediff -= rateShifts.get(i);
		return Math.exp(rates[i] - growth[i] * timediff);
	}

	private int getIntervalNr(double t) {
		// check which interval t + offset is in
		for (int i = 0; i < rateShifts.size()-1; i++)
			if (t < rateShifts.get(i+1))
				return i;
		// after the last interval, just keep using the last element
		return rateShifts.size()-1;
	}
	
	private int getLaterIntervalNr(double t, int startPoint) {
		// check which interval t + offset is in
		for (int i = startPoint; i < rateShifts.size()-1; i++)
			if (t < rateShifts.get(i+1))
				return i;
		// after the last interval, just keep using the last element
		return rateShifts.size()-1;
	}


	@Override
	public double getIntegral(double start, double finish) {
		if (start == finish)
			return 0.0;
		// get the interval "start" is in
		int first_int = getIntervalNr(start);
		// get the interval "finish" is in
		int last_int = getLaterIntervalNr(finish, first_int);
		
		double weighted = 0.0;
		double curr_time = start;	

		for (int i = first_int; i <= last_int; i++) {
			if (i > rateShifts.size()) {
				throw new IllegalArgumentException("rate shifts out of bounds");
			}

			double next_time = Math.min(getNextTime(i), finish);
			double r = growth[i];

			double rateShift = rateShifts.get(i);
			double timediff1 = curr_time - rateShift;
			double timediff2 = next_time - rateShift;

			
			if (r == 0.0) {
				weighted += (next_time - curr_time) * Math.exp(rates[i]);
			} else {
				weighted += invRatio[i] * ( Math.exp(-timediff2 * r) - Math.exp(-timediff1 * r));
			}

			curr_time = next_time;
		}
		return weighted;

	}
	
	
	@Override
	public double getIntensity(double t) {
		return getIntegral(0,t);
	}


	@Override
	public double getInverseIntensity(double x) {
		if (x == Double.POSITIVE_INFINITY) {
			return Double.POSITIVE_INFINITY;
        }
		
		int i = 0;
		double curr_time = 0;
		double integral = 0;
		
		do
		{
			double next_time = getNextTime(i);
			double r = growth[i];
	
//			double timediff1 = curr_time-rateShifts.get(i);
			double timediff2 = next_time-rateShifts.get(i);
	
			double old_diff = x - integral;
	
			if (r == 0.0) {
				integral += (next_time - curr_time) * Math.exp(rates[i]);
			} else {
				integral +=  invRatio[i] * (Math.exp(-timediff2 * r)-1.0);
			}
			
	
			double diff = x - integral;
	
			if (diff < 0 || i == rateShifts.size()) {
				
				if (r == 0.0) {
					return old_diff/Math.exp(rates[i]) + curr_time;
				} else {
					return Math.log((old_diff * -r)/ Math.exp(rates[i]) +1.0)/-r  + curr_time;
				}
			
			}
			
			curr_time = next_time;
			i++;
			
		}while(i<=rateShifts.size());
		
		
	
		return Double.POSITIVE_INFINITY;

	}


	private double getNextTime(int i) {
		if (i < rateShifts.size()-1)
			return rateShifts.get(i+1);
		else
			return Double.POSITIVE_INFINITY;
	}


	@Override
	public boolean requiresRecalculation() {
		recalculateNe();
		return super.requiresRecalculation();
	}

	@Override
	public void store() {
		growth_stored = new double[growth.length];
		System.arraycopy(growth, 0, growth_stored, 0, growth.length);
		stored_rates = new double[rates.length];
		System.arraycopy(rates, 0, stored_rates, 0, rates.length);
		invRatio_stored = new double[invRatio.length];
		System.arraycopy(invRatio, 0, invRatio_stored, 0, invRatio.length);
		super.store();
	}

	@Override
	public void restore() {
		System.arraycopy(growth_stored, 0, growth, 0, growth_stored.length);
		System.arraycopy(stored_rates, 0, rates, 0, stored_rates.length);
		System.arraycopy(invRatio_stored, 0, invRatio, 0, invRatio_stored.length);
		super.restore();
	}

	// computes the Ne's at the break points
	private void recalculateNe() {
		// logstandardize the active Predictor
		double[] logStandardPredictor = new double[rateShifts.size()];
		if (predictorIsActive.get()<predictors.size())  {	
			double mean = 0.0;
			for (int i = 0; i < predictors.size(); i++) {
			}
			for (int i = 0; i < independentAfterInput.get()+1; i++) {
				mean += predictors.get(predictorIsActive.get()).get(i);
			}
			mean /= (independentAfterInput.get()+1);
			for (int i = 0; i < independentAfterInput.get() + 1; i++) {
				logStandardPredictor[i] = predictors.get(predictorIsActive.get()).get(i) - mean;
			}
			double sd = 0.0;
			for (int i = 0; i < independentAfterInput.get() + 1; i++) {
				sd += Math.pow(logStandardPredictor[i], 2);
			}
			sd = Math.sqrt(sd / (independentAfterInput.get() + 1));
			for (int i = 0; i < independentAfterInput.get() + 1; i++) {
				logStandardPredictor[i] /= sd;
			}
		}
		

		
		growth = new double[rateShifts.size()];
		rates = new double[rateShifts.size()];
		double curr_time = 0.0;
		if (predictorIsActive.get()<predictors.size())  {			
			for (int i = 0; i < independentAfterInput.get()+1; i++) {
				rates[i] = effectSize.get(predictorIsActive.get())*
						logStandardPredictor[i] + NeToReassortment.get(i);
			}
			for (int i = independentAfterInput.get()+1; i < NeToReassortment.size(); i++) {
				rates[i] = NeToReassortment.get(i);
			}
			
		}else {
			for (int i = 0; i < NeToReassortment.size(); i++) {
				rates[i] = NeToReassortment.get(i);
			}
		}
		for (int i = 0; i < NeToReassortment.size()-1; i++) {
			growth[i] = (rates[i] - rates[i+1])
					/(rateShifts.get(i+1)-curr_time);
			curr_time = rateShifts.get(i+1);
			
			invRatio[i] = Math.exp(rates[i]) / -growth[i];
		}
		NesKnown = true;
	}

	@Override
	public void init(PrintStream out) {
		for (int i=0; i < rates.length; i++) {
			out.print("reassortmentRate." + i + "\t");
        }
		
	}

	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0; i < rates.length; i++) {
			out.print(stored_rates[i] + "\t");
		}
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
}