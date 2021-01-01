using System;
using System.Collections.Generic;

namespace CurveFit {
    public class CurveFitter : UserFunction {
        /** Constants for the built-in fit functions */
        public const int STRAIGHT_LINE = 0, POLY2 = 1, POLY3 = 2, POLY4 = 3,
                EXPONENTIAL = 4, POWER = 5, LOG = 6, RODBARD = 7, GAMMA_VARIATE = 8, LOG2 = 9,
                RODBARD2 = 10, EXP_WITH_OFFSET = 11, GAUSSIAN = 12, EXP_RECOVERY = 13,
                INV_RODBARD = 14, EXP_REGRESSION = 15, POWER_REGRESSION = 16,
                POLY5 = 17, POLY6 = 18, POLY7 = 19, POLY8 = 20,
                GAUSSIAN_NOOFFSET = 21, EXP_RECOVERY_NOOFFSET = 22, CHAPMAN = 23;

        /** Nicer sequence of the built-in function types */
        public static int[] sortedTypes = { STRAIGHT_LINE, POLY2, POLY3, POLY4,
            POLY5, POLY6, POLY7, POLY8,
            POWER, POWER_REGRESSION,
            EXPONENTIAL, EXP_REGRESSION, EXP_WITH_OFFSET,
            EXP_RECOVERY, EXP_RECOVERY_NOOFFSET,
            LOG, LOG2,
            GAUSSIAN, GAUSSIAN_NOOFFSET,
            RODBARD, RODBARD2, INV_RODBARD,
            GAMMA_VARIATE,CHAPMAN
    };

        /** Names of the built-in fit functions */
        public static string[] fitList = {"Straight Line","2nd Degree Polynomial",
            "3rd Degree Polynomial", "4th Degree Polynomial","Exponential","Power",
            "Log","Rodbard", "Gamma Variate", "y = a+b*ln(x-c)","Rodbard (NIH Image)",
            "Exponential with Offset","Gaussian", "Exponential Recovery",
            "Inverse Rodbard", "Exponential (linear regression)", "Power (linear regression)",
            "5th Degree Polynomial","6th Degree Polynomial","7th Degree Polynomial","8th Degree Polynomial",
            "Gaussian (no offset)", "Exponential Recovery (no offset)",
            "Chapman-Richards"
    }; // fList, doFit(), getNumParams() and makeInitialParamsAndVariations() must also be updated

        /** Equations of the built-in fit functions */
        public static string[] fList = {
            "y = a+bx","y = a+bx+cx^2",									//STRAIGHT_LINE,POLY2
            "y = a+bx+cx^2+dx^3","y = a+bx+cx^2+dx^3+ex^4",
            "y = a*exp(bx)","y = a*x^b", "y = a*ln(bx)",				//EXPONENTIAL,POWER,LOG
            "y = d+(a-d)/(1+(x/c)^b)", "y = b*(x-a)^c*exp(-(x-a)/d)",	//RODBARD,GAMMA_VARIATE
            "y = a+b*ln(x-c)", "x = d+(a-d)/(1+(y/c)^b) [y = c*((x-a)/(d-x))^(1/b)]",  //LOG2,RODBARD2
            "y = a*exp(-bx) + c", "y = a + (b-a)*exp(-(x-c)*(x-c)/(2*d*d))", //EXP_WITH_OFFSET,GAUSSIAN
            "y = a*(1-exp(-b*x)) + c", "y = c*((x-a)/(d-x))^(1/b)",		//EXP_RECOVERY, INV_RODBARD
            "y = a*exp(bx)", "y = a*x^b",								//EXP_REGRESSION, POWER_REGRESSION
            "y = a+bx+cx^2+dx^3+ex^4+fx^5", "y = a+bx+cx^2+dx^3+ex^4+fx^5+gx^6",
            "y = a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7", "y = a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7+ix^8",
            "y = a*exp(-(x-b)*(x-b)/(2*c*c)))",							//GAUSSIAN_NOOFFSET
            "y = a*(1-exp(-b*x))",                                      //EXP_RECOVERY_NOOFFSET
            "y = a*(1-exp(-b*x))^c"									    //CHAPMAN
    };

        /** @deprecated now in the io.binroot.regression.Minimizer class (since ImageJ 1.46f).
         *	(probably of not much value for anyone anyhow?) */
        public static int IterFactor = 500;

        private const int CUSTOM = 100;   // user function defined in macro or plugin
        private const int GAUSSIAN_INTERNAL = 101;   // Gaussian with separate offset & multiplier
        private const int RODBARD_INTERNAL = 102;    // Rodbard with separate offset & multiplier

        private int fitType = -1;       // Number of curve type to fit
        private double[] xData, yData;  // x,y data to fit
        private double[] xDataSave, yDataSave;  //saved original data after fitting modified data
        private int numPoints;          // number of data points in actual fit
        private double ySign = 0;       // remember sign of y data for power-law fit via regression
        private double sumY = Double.NaN, sumY2 = Double.NaN;  // sum(y), sum(y^2) of the data used for fitting
        private int numParams;          // number of parameters
        private double[] initialParams; // user specified or estimated initial parameters
        private double[] initialParamVariations; // estimate of range of parameters
        private double[] minimizerInitialParams;          // modified initialParams of the minimizer
        private double[] minimizerInitialParamVariations; // modified initialParamVariations of the minimizer
        private double maxRelError = 1e-10;// maximum relative deviation of sum of residuals^2 from minimum
        private long time;              // elapsed time in ms
        private int customParamCount;   // number of parameters of user-supplied function (macro or plugin)
        private String customFormula;   // equation defined in macro or text
        private UserFunction userFunction; // plugin with custom fit function
        private int macroStartProgramCounter;    // equation in macro starts here
        private int numRegressionParams;// number of parameters that can be calculated by linear regression
        private int offsetParam = -1;   // parameter number: function has this parameter as offset
        private int factorParam = -1;   // index of parameter that is slope or multiplicative factor of the function
        private bool hasSlopeParam;  // whether the 'factorParam' is the slope of the function
        private double[] finalParams;   // parameters obtained by fit
        private bool linearRegressionUsed;   // whether linear regression alone was used instead of the minimizer
        private bool restrictPower;  // power via linear regression fit: (0,0) requires positive power
        private Minimizer minimizer = new Minimizer();
        private int minimizerStatus = Minimizer.INITIALIZATION_FAILURE; // status of the minimizer after minimizing
        private String errorString;     // in case of error before invoking the minimizer
        private static String[] sortedFitList; // names like fitList, but in more logical sequence
        private static Dictionary<string, int> namesTable; // converts fitList String into number

        /** Construct a new io.binroot.regression.CurveFitter. */
        public CurveFitter(double[] xData, double[] yData) {
            this.xData = xData;
            this.yData = yData;
            numPoints = xData.Length;
        }

        /** Perform curve fitting with one of the built-in functions
         *			doFit(fitType) does the fit quietly
         *	Use getStatus() and/or getStatusString() to see whether fitting was (probably) successful and
         *	getParams() to access the result.
         */
        public void doFit(int fitType) {
            doFit(fitType, false);
        }

        /** Perform curve fitting with one of the built-in functions
         *			doFit(fitType, true) pops up a dialog allowing the user to set the initial
         *						fit parameters and various numbers controlling the io.binroot.regression.Minimizer
         *	Use getStatus() and/or getStatusString() to see whether fitting was (probably) successful and
         *	getParams() to access the result.
         */
        public void doFit(int fitType, bool showSettings) {
            if (!(fitType >= STRAIGHT_LINE && fitType < fitList.Length || fitType == CUSTOM))
                throw new ArgumentOutOfRangeException("Invalid fit type");
            if (fitType == CUSTOM && userFunction == null)
                throw new ArgumentOutOfRangeException("No custom formula!");
            this.fitType = fitType;
            if (isModifiedFitType(fitType))         // these fits don't use the original data and a different fit type (this.fitType)
                if (!prepareModifiedFitType(fitType)) return;
            numParams = getNumParams();
            if (fitType != CUSTOM)
                getOffsetAndFactorParams();
            //io.binroot.regression.IJ.log("special params1: off="+offsetParam+(hasSlopeParam ? " slo=" : " fac=")+factorParam+" numPar="+numParams+" numRegressPar="+numRegressionParams);
            calculateSumYandY2();                   // sumY, sumY2 needed for regression, abs Error; R, goodness of modified fit functions
            long startTime = Java.SystemCurrentTimeMillis();
            if (this.fitType == STRAIGHT_LINE) {    // no minimizer needed
                finalParams = new double[] { 0, 0, 0 }; // (does not work with 0 initial slope)
                doRegression(finalParams);          // fit by regression; save sum(Residuals^2) as last array element
                linearRegressionUsed = true;
            } else {                                // use simplex minimizer
                minimizer.setFunction(this, numParams - numRegressionParams);
                minimizer.setExtraArrayElements(numRegressionParams);  // reserve space for extra parameters if minimizer has fewer paramters
                if (!makeInitialParamsAndVariations(fitType))       // also includes some data checking
                    return;                         // initialization failure
                if (numRegressionParams > 0)
                    modifyInitialParamsAndVariations();
                else {
                    minimizerInitialParams = initialParams;
                    minimizerInitialParamVariations = initialParamVariations;
                }
                startTime = Java.SystemCurrentTimeMillis();
                // The maximum absolute error of the fit must be specified in case the
                // fit function fits perfectly, i.e. the sume of residuals approaches 0.
                // In such a case, the maximum relative error is meaningless and the
                // minimizer would run until it reaches the maximum iteration count.
                double maxAbsError = Math.Min(1e-6, maxRelError) * Math.Sqrt(sumY2);
                minimizer.setMaxError(maxRelError, maxAbsError);
                //{String s="initVariations:";for(int ii=0;ii<numParams;ii++)s+=" ["+ii+"]:"+io.binroot.regression.IJ.d2s(initialParamVariations[ii],5,9);io.binroot.regression.IJ.log(s);}
                //{String s="minInitVariations:";for(int ii=0;ii<numParams;ii++)s+=" ["+ii+"]:"+io.binroot.regression.IJ.d2s(minimizerInitialParamVariations[ii],5,9);io.binroot.regression.IJ.log(s);}
                //{String s="minInitPars:";for(int ii=0;ii<numParams;ii++)s+=" ["+ii+"]:"+io.binroot.regression.IJ.d2s(minimizerInitialParams[ii],5,9);io.binroot.regression.IJ.log(s);}
                // m i n i m i z a t i o n	of squared residuals
                minimizerStatus = minimizer.minimize(minimizerInitialParams, minimizerInitialParamVariations);
                finalParams = minimizer.getParams();
                if (numRegressionParams > 0)
                    minimizerParamsToFullParams(finalParams, false);
            }
            if (isModifiedFitType(fitType))         //params1 of actual fit to user params1
                postProcessModifiedFitType(fitType);

            switch (fitType) {                      //postprocessing for nicer output:
                case GAUSSIAN:                      //Gaussians: width (std deviation) should be >0
                    finalParams[3] = Math.Abs(finalParams[3]); break;
                case GAUSSIAN_NOOFFSET:
                    finalParams[2] = Math.Abs(finalParams[2]); break;
            }
            time = Java.SystemCurrentTimeMillis() - startTime;
        }

        /** Sets the initial parameters, which override the default initial parameters. */
        public void setInitialParameters(double[] initialParams) {
            this.initialParams = initialParams;
        }

        /** Returns a reference to the io.binroot.regression.Minimizer used, for accessing io.binroot.regression.Minimizer methods directly.
         *	Note that no io.binroot.regression.Minimizer is used if fitType is any of STRAIGHT_LINE, EXP_REGRESSION,
         *	and POWER_REGRESSION. */
        public Minimizer getMinimizer() {
            return minimizer;
        }

        /** For improved fitting performance when using a custom fit formula, one may
         *	specify parameters that can be calculated directly by linear regression.
         *	For values not used, set the index to -1
         *
         * @param offsetParam  Index of a parameter that is a pure offset:
         *					   E.g. '0' if	f(p0, p1, p2...) = p0 + function(p1, p2, ...).
         * @param multiplyParam	 Index of a parameter that is purely multiplicative.
         *					   E.g. multiplyParams=1 if f(p0, p1, p2, p3...) can be expressed as
         *					   p1*func(p0, p2, p3, ...) or p0 +p1*func(p0, p2, p3, ...) with '0' being
         *					   the offsetparam.
         * @param slopeParam   Index of a parameter that is multiplied with x and then summed to the function.
         *					   E.g. '1' for f(p0, p1, p2, p3...) = p1*x + func(p0, p2, p3, ...)
         *					   Only one, multiplyParam and slopeParam can be used (ie.e, the other
         *					   should be set to -1)
         */
        public void setOffsetMultiplySlopeParams(int offsetParam, int multiplyParam, int slopeParam) {
            this.offsetParam = offsetParam;
            hasSlopeParam = slopeParam >= 0;
            factorParam = hasSlopeParam ? slopeParam : multiplyParam;
            numRegressionParams = 0;
            if (offsetParam >= 0) numRegressionParams++;
            if (factorParam >= 0) numRegressionParams++;
        }

        /** Get number of parameters for current fit formula
         *	Do not use before 'doFit', because the fit function would be undefined.	 */
        public int getNumParams() {
            switch (fitType) {
                case STRAIGHT_LINE: return 2;
                case POLY2: return 3;
                case POLY3: return 4;
                case POLY4: return 5;
                case POLY5: return 6;
                case POLY6: return 7;
                case POLY7: return 8;
                case POLY8: return 9;
                case EXPONENTIAL: case EXP_REGRESSION: return 2;
                case POWER: case POWER_REGRESSION: return 2;
                case EXP_RECOVERY_NOOFFSET: return 2;
                case LOG: return 2;
                case LOG2: return 3;
                case GAUSSIAN_NOOFFSET: return 3;
                case EXP_RECOVERY: return 3;
                case CHAPMAN: return 3;
                case EXP_WITH_OFFSET: return 3;
                case RODBARD: case RODBARD2: case INV_RODBARD: case RODBARD_INTERNAL: return 4;
                case GAMMA_VARIATE: return 4;
                case GAUSSIAN: case GAUSSIAN_INTERNAL: return 4;
                case CUSTOM: return customParamCount;
            }
            return 0;
        }

        /** Returns the formula value for parameters 'p' at 'x'.
         *	Do not use before 'doFit', because the fit function would be undefined.	 */
        public double f(double[] p, double x) {
            if (fitType != CUSTOM)
                return f(fitType, p, x);
            else {
                return userFunction.userFunction2(p, x);
            }
        }

        /** Returns value of built-in 'fitType' formula value for parameters "p" at "x" */
        public static double f(int fitType, double[] p, double x) {
            switch (fitType) {
                case STRAIGHT_LINE:
                    return p[0] + x * p[1];
                case POLY2:
                    return p[0] + x * (p[1] + x * p[2]);
                case POLY3:
                    return p[0] + x * (p[1] + x * (p[2] + x * p[3]));
                case POLY4:
                    return p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * p[4])));
                case POLY5:
                    return p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * p[5]))));
                case POLY6:
                    return p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * (p[5] + x * p[6])))));
                case POLY7:
                    return p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * (p[5] + x * (p[6] + x * p[7]))))));
                case POLY8:
                    return p[0] + x * (p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * (p[5] + x * (p[6] + x * (p[7] + x * p[8])))))));
                case EXPONENTIAL:
                case EXP_REGRESSION:
                    return p[0] * Math.Exp(p[1] * x);
                case EXP_WITH_OFFSET:
                    return p[0] * Math.Exp(-p[1] * x) + p[2];
                case EXP_RECOVERY:
                    return p[0] * (1 - Math.Exp(-p[1] * x)) + p[2];
                case EXP_RECOVERY_NOOFFSET:
                    return p[0] * (1 - Math.Exp(-p[1] * x));
                case CHAPMAN:                               // a*(1-exp(-b*x))^c
                    double value = p[0] * (Math.Pow((1 - (Math.Exp(-p[1] * x))), p[2]));
                    //	Log.e("test", "values = " + value);
                    return value;
                case GAUSSIAN:
                    return p[0] + (p[1] - p[0]) * Math.Exp(-(x - p[2]) * (x - p[2]) / (2.0 * p[3] * p[3]));
                case GAUSSIAN_INTERNAL:                     // replaces GAUSSIAN for the fitting process
                    return p[0] + p[1] * Math.Exp(-(x - p[2]) * (x - p[2]) / (2.0 * p[3] * p[3]));
                case GAUSSIAN_NOOFFSET:
                    return p[0] * Math.Exp(-(x - p[1]) * (x - p[1]) / (2.0 * p[2] * p[2]));
                case POWER:                                 // ax^b
                case POWER_REGRESSION:
                    return p[0] * Math.Pow(x, p[1]);
                case LOG:
                    if (x == 0.0)
                        return -1000 * p[0];
                    return p[0] * Math.Log(p[1] * x);
                case RODBARD: {                             // d+(a-d)/(1+(x/c)^b)
                        double ex = Math.Pow(x / p[2], p[1]); //(x/c)^b
                        return p[3] + (p[0] - p[3]) / (1.0 + ex);
                    }
                case RODBARD_INTERNAL: {                    // d+a/(1+(x/c)^b) , replaces RODBARD of the fitting process
                        double ex = Math.Pow(x / p[2], p[1]); //(x/c)^b
                        return p[3] + p[0] / (1.0 + ex);
                    }
                case GAMMA_VARIATE:                         // b*(x-a)^c*exp(-(x-a)/d)
                    if (p[0] >= x) return 0.0;
                    if (p[1] <= 0) return Double.NaN;
                    if (p[2] <= 0) return Double.NaN;
                    if (p[3] <= 0) return Double.NaN;

                    double pw = Math.Pow((x - p[0]), p[2]);
                    double e = Math.Exp((-(x - p[0])) / p[3]);
                    return p[1] * pw * e;
                case LOG2:
                    double tmp = x - p[2];
                    if (tmp <= 0)
                        return double.NaN;
                    return p[0] + p[1] * Math.Log(tmp);
                case INV_RODBARD:       // c*((x-a)/(d-x))^(1/b), the inverse Rodbard function
                case RODBARD2:          // also used after the 'Rodbard NIH Image' fit
                    double y;
                    if (p[3] - x < 2 * double.MinValue || x < p[0]) // avoid x>=d (singularity) and x<a (negative exponent)
                        y = fitType == INV_RODBARD ? Double.NaN : 0.0;
                    else {
                        y = (x - p[0]) / (p[3] - x); //(x-a)/(d-x) = ( (a-d)/(x-d) -1 )
                        y = Math.Pow(y, 1.0 / p[1]);  //y=y**(1/b)
                        y = y * p[2];
                    }
                    return y;
                default:
                    return 0.0;
            }
        }

        /** Get the result of fitting, i.e. the set of parameter values for the best fit.
         *	Note that the array returned may have more elements than numParams; ignore the rest.
         *	May return an array with only NaN values if the minimizer could not start properly,
         *	i.e., if getStatus() returns io.binroot.regression.Minimizer.INITILIZATION_FAILURE.
         *	See io.binroot.regression.Minimizer.getParams() for details.
         */
        public double[] getParams() {
            return finalParams == null ? minimizer.getParams() : finalParams; //if we have no result, take all_NaN result from the io.binroot.regression.Minimizer
        }

        /** Returns residuals array, i.e., differences between data and curve.
         *	The residuals are with respect to the real data, also for fit types where the data are
         *	modified before fitting (power&exp fit by linear regression, 'Rodbard NIH Image' ).
         *	This is in contrast to sum of squared residuals, which is for the fit that was actually done.
         */
        public double[] getResiduals() {
            double[] params1 = getParams();
            double[] residuals = new double[xData.Length];
            for (int i = 0; i < xData.Length; i++)
                residuals[i] = yData[i] - f(params1, xData[i]);
            return residuals;
        }

        /* Get the sum of the residuals (may be NaN if the minimizer could not start properly
         *	i.e., if getStatus() returns io.binroot.regression.Minimizer.INITILIZATION_FAILURE).
         */
        public double getSumResidualsSqr() {
            return getParams()[numParams];  // value is stored as last element by the minimizer
        }

        /** Returns the standard deviation of the residuals.
         *	Here, the standard deviation is defined here as the root-mean-square of the residuals
         *	times sqrt(n/(n-1)); where n is the number of points.
         */
        public double getSD() {
            double[] residuals = getResiduals();
            int n = residuals.Length;
            double sum = 0.0, sum2 = 0.0;
            for (int i = 0; i < n; i++) {
                sum += residuals[i];
                sum2 += residuals[i] * residuals[i];
            }
            double stdDev = (sum2 - sum * sum / n); //sum of squared residuals
            return Math.Sqrt(stdDev / (n - 1.0));
        }

        /** Returns R^2, where 1.0 is best.
         <pre>
         r^2 = 1 - SSE/SSD

         where:	 SSE = sum of the squared errors
         SSD = sum of the squared deviations about the mean.
         </pre>
         *	For power, exp by linear regression and 'Rodbard NIH Image', this is calculated for the
         *	fit actually done, not for the residuals of the original data.
         */
        public double getRSquared() {
            if (double.IsNaN(sumY)) calculateSumYandY2();
            double sumMeanDiffSqr = sumY2 - sumY * sumY / numPoints;
            double rSquared = 0.0;
            if (sumMeanDiffSqr > 0.0)
                rSquared = 1.0 - getSumResidualsSqr() / sumMeanDiffSqr;
            return rSquared;
        }

        /** Get a measure of "goodness of fit" where 1.0 is best.
         *	Approaches R^2 if the number of points is much larger than the number of fit parameters.
         *	For power, exp by linear regression and 'Rodbard NIH Image', this is calculated for the
         *	fit actually done, not for the residuals of the original data.
         */
        public double getFitGoodness() {
            if (double.IsNaN(sumY)) calculateSumYandY2();
            double sumMeanDiffSqr = sumY2 - sumY * sumY / numPoints;
            double fitGoodness = 0.0;
            int degreesOfFreedom = numPoints - getNumParams();
            if (sumMeanDiffSqr > 0.0 && degreesOfFreedom > 0)
                fitGoodness = 1.0 - (getSumResidualsSqr() / sumMeanDiffSqr) * numPoints / (double)degreesOfFreedom;

            return fitGoodness;
        }

        public int getStatus() {
            return linearRegressionUsed ? Minimizer.SUCCESS : minimizerStatus;
        }

        /** Get a short text with a short description of the status. Should be preferred over
         *	io.binroot.regression.Minimizer.STATUS_STRING[getMinimizer().getStatus()] because getStatusString()
         *	better explains the problem in some cases of initialization failure (data not
         *	compatible with the fit function chosen) */
        public String getStatusString() {
            return errorString != null ? errorString : Minimizer.STATUS_STRING[getStatus()];
        }

        /** Get a string with detailed description of the curve fitting results (several lines,
         *	including the fit parameters).
         */
        public String getResultString() {
            String resultS = "\nFormula: " + getFormula() +
                    "\nStatus: " + getStatusString();
            if (!linearRegressionUsed) resultS += "\nNumber of completed minimizations: " + minimizer.getCompletedMinimizations();
            resultS += "\nNumber of iterations: " + getIterations();
            if (!linearRegressionUsed) resultS += " (max: " + minimizer.getMaxIterations() + ")";
            resultS += "\nTime: " + time + " ms" +
                    "\nSum of residuals squared: " + IJ.d2s(getSumResidualsSqr(), 5, 9) +
                    "\nStandard deviation: " + IJ.d2s(getSD(), 5, 9) +
                    "\nR^2: " + IJ.d2s(getRSquared(), 5) +
                    "\nParameters:";
            char pChar = 'a';
            double[] pVal = getParams();
            for (int i = 0; i < numParams; i++) {
                resultS += "\n	" + pChar + " = " + IJ.d2s(pVal[i], 5, 9);
                pChar++;
            }
            return resultS;
        }

        /** Set maximum number of simplex restarts to do. See io.binroot.regression.Minimizer.setMaxRestarts for details. */
        public void setRestarts(int maxRestarts) {
            minimizer.setMaxRestarts(maxRestarts);
        }

        /** Set the maximum error. by which the sum of residuals may deviate from the true value
         *	(relative w.r.t. full sum of rediduals). Possible range: 0.1 ... 10^-16 */
        public void setMaxError(double maxRelError) {
            if (double.IsNaN(maxRelError)) return;
            if (maxRelError > 0.1) maxRelError = 0.1;
            if (maxRelError < 1e-16) maxRelError = 1e-16;   // can't be less than numerical accuracy
            this.maxRelError = maxRelError;
        }

        /** Create output on the number of iterations in the ImageJ Status line, e.g.
         *  "<ijStatusString> 50 (max 750); ESC to stop"
         *  @param ijStatusString Displayed in the beginning of the status message. No display if null.
         *  E.g. "Curve Fit: Iteration "
         *  @param checkEscape When true, the io.binroot.regression.Minimizer stops if escape is pressed and the status
         *  becomes ABORTED. Note that checking for ESC does not work in the Event Queue thread. */
        public void setStatusAndEsc(String ijStatusString, bool checkEscape) {
            minimizer.setStatusAndEsc(ijStatusString, checkEscape);
        }

        /** Get number of iterations performed. Returns 1 in case the fit was done by linear regression only. */
        public int getIterations() {
            return linearRegressionUsed ? 1 : minimizer.getIterations();
        }

        /** Get maximum number of iterations allowed (sum of iteration count for all restarts) */
        public int getMaxIterations() {
            return minimizer.getMaxIterations();
        }

        /** Set maximum number of iterations allowed (sum of iteration count for all restarts) */
        public void setMaxIterations(int maxIter) {
            minimizer.setMaxIterations(maxIter);
        }

        /** Get maximum number of simplex restarts to do. See io.binroot.regression.Minimizer.setMaxRestarts for details. */
        public int getRestarts() {
            return minimizer.getMaxRestarts();
        }

        /** Returns the status of the io.binroot.regression.Minimizer after doFit.  io.binroot.regression.Minimizer.SUCCESS indicates a
         *	successful completion. In case of io.binroot.regression.Minimizer.INITIALIZATION_FAILURE, fitting could
         *	not be performed because the data and/or initial parameters are not compatible
         *	with the function value.  In that case, getStatusString may explain the problem.
         *	For further status codes indicating problems during fitting, see the status codes
         *	of the Minimzer class. */

        /** returns the array with the x data */
        public double[] getXPoints() {
            return xData;
        }

        /** returns the array with the y data */
        public double[] getYPoints() {
            return yData;
        }

        /** returns the code of the fit type of the fit performed */
        public int getFit() {
            return fitType;
        }

        /** returns the name of the fit function of the fit performed */
        public String getName() {
            if (fitType == CUSTOM)
                return "User-defined";
            if (fitType == GAUSSIAN_INTERNAL)
                fitType = GAUSSIAN;
            else if (fitType == RODBARD_INTERNAL)
                fitType = RODBARD;
            return fitList[fitType];
        }

        /** returns a String with the formula of the fit function used */
        public String getFormula() {
            if (fitType == CUSTOM)
                return customFormula;
            else
                return fList[fitType];
        }

        /** Returns an array of fit names with nicer sorting */
        public static String[] getSortedFitList() {
            if (sortedFitList == null) {
                String[] l = new String[fitList.Length];
                for (int i = 0; i < fitList.Length; i++)
                    sortedFitList[i] = fitList[sortedTypes[i]];
                sortedFitList = l;
            }
            return sortedFitList;
        }

        /** Returns the code for a fit with given name as defined in fitList, or -1 if not found */
        public static int getFitCode(String fitName) {
            if (namesTable == null) {
                var h = new Dictionary<String, int>();
                for (int i = 0; i < fitList.Length; i++)
                    h.Add(fitList[i], i);
                namesTable = h;
            }
            return namesTable[fitName];
        }

        /** This function is called by the io.binroot.regression.Minimizer and calculates the sum of squared
         *	residuals for given parameters.
         *	To improve the efficiency, simple linear dependencies are solved directly
         *	by linear regression; in that case the corresponding parameters are modified.
         *	This effectively reduces the number of free parameters by one or two and thereby
         *	significantly improves the performance of minimization.
         */
        public double userFunction2(double[] params1, double dummy) {
            double sumResidualsSqr = 0.0;
            if (numRegressionParams == 0) {     // simply calculate sum of residuals
                for (int i = 0; i < numPoints; i++) {
                    double fValue = f(params1, xData[i]);
                    sumResidualsSqr += sqr(fValue - yData[i]);
                }
                //io.binroot.regression.IJ.log(io.binroot.regression.IJ.d2s(params1[0],3,5)+","+io.binroot.regression.IJ.d2s(params1[1],3,5)+": r="+io.binroot.regression.IJ.d2s(sumResidualsSqr,3,5)+Thread.currentThread().getName() );

            } else {    // handle simple linear dependencies by linear regression:
                        //if(getIterations()<1){String s="minimizerPar:";for(int ii=0;ii<=numParams;ii++)s+=" ["+ii+"]:"+io.binroot.regression.IJ.d2s(params1[ii],5,9);io.binroot.regression.IJ.log(s);}
                minimizerParamsToFullParams(params1, true);
                doRegression(params1);
                sumResidualsSqr = fullParamsToMinimizerParams(params1);
            }
            return sumResidualsSqr;
        }

        /** For fits where linear regression is used to reduce the number of parameters handled
         *	by the io.binroot.regression.Minimizer, convert io.binroot.regression.Minimizer parameters to the complete set of parameters.
         *	When not for calculating regression, we use the sum of squared residuals,
         *	offset and multiplication factor stored in the extra array elements:
         *	The io.binroot.regression.Minimizer stores the sum of squared residuals directly behind its last parameter.
         *	The next element is the value of the offsetParam (if any).
         *	The element is the value of the factorParam (if any).
         */
        private void minimizerParamsToFullParams(double[] params1, bool forRegression) {
            bool shouldTransformToSmallerParams = false;
            double offset = 0;
            double factor = hasSlopeParam ? 0 : 1; //for regression, calculate function value without x*slope, but multiplied with 1
            double sumResidualsSqr = 0;
            if (!forRegression) {               // recover regression-calculated parameters from extra array elements
                int i = params1.Length - 1;
                if (factorParam >= 0)
                    factor = params1[i--];
                if (offsetParam >= 0)
                    offset = params1[i];
                sumResidualsSqr = params1[numParams - numRegressionParams]; // sum of squared residuals has been calculated already
                params1[numParams] = sumResidualsSqr;                       // ... and is now stored in its new (proper) place
            }
            // move array elements to position in array with full number of parameters
            for (int i = numParams - 1, iM = numParams - numRegressionParams - 1; i >= 0; i--) {
                if (i == offsetParam)
                    params1[i] = offset;
                else if (i == factorParam)
                    params1[i] = factor;
                else
                    params1[i] = params1[iM--];
            }
            params1[numParams] = sumResidualsSqr;
        }

        /** Determine sum of squared residuals with linear regression.
         *	The sum of squared residuals is written to the array element with index 'numParams',
         *	the offset and factor params1 (if any) are written to their proper positions in the
         *	params1 array */
        private void doRegression(double[] params1) {
            double sumX = 0, sumX2 = 0, sumXY = 0; //sums for regression; here 'x' are function values
            double sumY = 0, sumY2 = 0;         //only calculated for 'slope', otherwise we use the values calculated already
            for (int i = 0; i < numPoints; i++) {
                double fValue = fitType == STRAIGHT_LINE ? 0 : f(params1, xData[i]);  // function value
                if (double.IsNaN(fValue)) { //check for NaN now; later we need NaN checking for division-by-zero check.
                    params1[numParams] = Double.NaN;
                    return;                 //sum of squared residuals is NaN if any value is NaN
                }
                //if(getIterations()==0)io.binroot.regression.IJ.log(xData[i]+"\t"+yData[i]+"\t"+fValue); //x,y,function
                if (hasSlopeParam) {        // fit y = offset + slope*x + function(of other params1)
                    double x = xData[i];
                    double y = yData[i] - fValue;
                    sumX += x;
                    sumX2 += x * x;
                    sumXY += x * y;
                    sumY2 += y * y;
                    sumY += y;
                } else {                    // fit y = offset + factor * function(of other params1)
                    double x = fValue;
                    double y = yData[i];
                    sumX += fValue;
                    sumX2 += fValue * fValue;
                    sumXY += fValue * yData[i];
                }
            }
            if (!hasSlopeParam) {
                sumY = this.sumY;
                sumY2 = this.sumY2;
            }
            double factor = 0; // factor or slope
            double sumResidualsSqr = 0;
            if (offsetParam < 0) {          // only 'factor' parameter, no offset:
                factor = sumXY / sumX2;
                if (double.IsNaN(factor) || double.IsInfinity(factor))
                    factor = 0;             // all 'x' values are 0, any factor (slope) will fit
                sumResidualsSqr = sumY2 + factor * (factor * sumX2 - 2 * sumXY);
                if (sumResidualsSqr < 2e-15 * sumY2)
                    sumResidualsSqr = 2e-15 * sumY2;
            } else {                        // full linear regression or offset only. Slope is named 'factor' here
                if (factorParam >= 0) {
                    factor = (sumXY - sumX * sumY / numPoints) / (sumX2 - sumX * sumX / numPoints);
                    if (restrictPower & factor <= 0)    // power-law fit with (0,0) point: power must be >0
                        factor = 1e-100;
                    else if (double.IsNaN(factor) || double.IsInfinity(factor))
                        factor = 0;         // all 'x' values are equal, any factor (slope) will fit
                }
                double offset = (sumY - factor * sumX) / numPoints;
                params1[offsetParam] = offset;
                sumResidualsSqr = sqr(factor) * sumX2 + numPoints * sqr(offset) + sumY2 +
                        2 * factor * offset * sumX - 2 * factor * sumXY - 2 * offset * sumY;
                // check for accuracy problem: large difference of small numbers?
                // Don't report unrealistic or even negative values, otherwise minimization could lead
                // into parameters where we have a numeric problem
                if (sumResidualsSqr < 2e-15 * (sqr(factor) * sumX2 + numPoints * sqr(offset) + sumY2))
                    sumResidualsSqr = 2e-15 * (sqr(factor) * sumX2 + numPoints * sqr(offset) + sumY2);
                //if(){io.binroot.regression.IJ.log("sumX="+sumX+" sumX2="+sumX2+" sumXY="+sumXY+" factor="+factor+" offset=="+offset);}
            }
            params1[numParams] = sumResidualsSqr;
            if (factorParam >= 0)
                params1[factorParam] = factor;
        }
        /** Convert full set of parameters to minimizer parameters and returns the sum of residuals squared.
         *	The last array elements, not used by the minimizer, are the offset and factor parameters (if any)
         */
        private double fullParamsToMinimizerParams(double[] params1) {
            double offset = offsetParam >= 0 ? params1[offsetParam] : 0;
            double factor = factorParam >= 0 ? params1[factorParam] : 0;
            double sumResidualsSqr = params1[numParams];

            for (int x = 0, iNew = 0; x < numParams; x++) {       // leave only the parameters for the minimizer in the beginning of the array
                if (x != factorParam && x != offsetParam) {
                    params1[iNew++] = params1[x];
                }
            }
            int i = params1.Length - 1;
            if (factorParam >= 0)
                params1[i--] = factor;
            if (offsetParam >= 0)
                params1[i--] = offset;
            params1[i--] = sumResidualsSqr;
            return sumResidualsSqr;
        }

        /** In case one or two parameters are calculated by regression and not by the minimizer:
         *	Make modified initialParams and initialParamVariations for the io.binroot.regression.Minimizer
         */
        private void modifyInitialParamsAndVariations() {
            minimizerInitialParams = (double[])initialParams.Clone();
            minimizerInitialParamVariations = (double[])initialParamVariations.Clone();
            if (numRegressionParams > 0) // convert to shorter arrays with only the parameters used by the minimizer
                for (int i = 0, iNew = 0; i < numParams; i++)
                    if (i != factorParam && i != offsetParam) {
                        minimizerInitialParams[iNew] = minimizerInitialParams[i];
                        minimizerInitialParamVariations[iNew] = minimizerInitialParamVariations[i];
                        iNew++;
                    }
        }

        /** Estimate values for initial parameters and their typical range for the built-in
         *	function types.	 For fits with modifications for regression, 'fitType' is still
         *	the type of the original (unmodified) fit.
         *	Also checks for x values being non-negative for fit types that require this,
         *	and returns false if the data cannot be fitted because of negative x.
         *	In such a case, 'errorString' contains a message about the problem. */
        private bool makeInitialParamsAndVariations(int fitType) {
            bool hasInitialParams = initialParams != null;
            bool hasInitialParamVariations = initialParamVariations != null;
            if (!hasInitialParams) {
                initialParams = new double[numParams];
                if (fitType == CUSTOM) {
                    for (int i = 0; i < numParams; i++)
                        initialParams[i] = 1.0;
                }
            }
            if (fitType == CUSTOM)
                return true; // no way to guess initial parameters or initialParamVariations
            if (!hasInitialParamVariations)
                initialParamVariations = new double[numParams];

            // Calculate some things that might be useful for predicting parameters
            double firstx = xData[0];
            double firsty = yData[0];
            double lastx = xData[numPoints - 1];
            double lasty = yData[numPoints - 1];
            double xMin = firstx, xMax = firstx;    //with this initialization, the loop starts at 1, not 0
            double yMin = firsty, yMax = firsty;
            double xMean = firstx, yMean = firsty;
            double xOfMax = firstx;
            for (int i = 1; i < numPoints; i++) {
                xMean += xData[i];
                yMean += yData[i];
                if (xData[i] > xMax) xMax = xData[i];
                if (xData[i] < xMin) xMin = xData[i];
                if (yData[i] > yMax) { yMax = yData[i]; xOfMax = xData[i]; }
                if (yData[i] < yMin) yMin = yData[i];
            }
            xMean /= numPoints;
            yMean /= numPoints;

            double slope = (lasty - firsty) / (lastx - firstx);
            if (double.IsNaN(slope) || double.IsInfinity(slope)) slope = 0;
            double yIntercept = yMean - slope * xMean;

            //We cannot fit the following cases because we would always get NaN function values
            if (xMin < 0 && (fitType == POWER || fitType == CHAPMAN)) {
                errorString = "Cannot fit " + fitList[fitType] + " when x<0";
                return false;
            } else if (xMin < 0 && xMax > 0 && fitType == RODBARD) {
                errorString = "Cannot fit " + fitList[fitType] + " to mixture of x<0 and x>0";
                return false;
            } else if (xMin <= 0 && fitType == LOG) {
                errorString = "Cannot fit " + fitList[fitType] + " when x<=0";
                return false;
            }

            if (!hasInitialParams) {
                switch (fitType) {
                    //case POLY2: case POLY3: case POLY4: case POLY5: case POLY6: case POLY7: case POLY8:
                    // offset&slope calculated via regression; leave the others at 0

                    // also for the other cases, some initial parameters are unused; only to show them with 'showSettings'
                    case EXPONENTIAL:           // a*exp(bx)   assuming change by factor of e between xMin & xMax
                        initialParams[1] = 1.0 / (xMax - xMin + 1e-100) * Math.Sign(yMean) * Math.Sign(slope);
                        initialParams[0] = yMean / Math.Exp(initialParams[1] * xMean); //don't care, done by regression
                        break;
                    case EXP_WITH_OFFSET:       // a*exp(-bx) + c	assuming b>0, change by factor of e between xMin & xMax
                    case EXP_RECOVERY:          // a*(1-exp(-bx)) + c
                        initialParams[1] = 1.0 / (xMax - xMin + 1e-100);
                        initialParams[0] = (yMax - yMin) / Math.Exp(initialParams[1] * xMean) * Math.Sign(slope) *
                                fitType == EXP_RECOVERY ? 1 : -1; // don't care, we will do this via regression
                        initialParams[2] = 0.5 * yMean;         // don't care, we will do this via regression
                        break;
                    case EXP_RECOVERY_NOOFFSET: // a*(1-exp(-bx))
                        initialParams[1] = 1.0 / (xMax - xMin + 1e-100) * Math.Sign(yMean) * Math.Sign(slope);
                        initialParams[0] = yMean / Math.Exp(initialParams[1] * xMean); //don't care, done by regression
                        break;
                    case POWER:                 // ax^b, assume linear for the beginning
                        initialParams[0] = yMean / (Math.Abs(xMean + 1e-100));  // don't care, we will do this via regression
                        initialParams[1] = 1.0;
                        break;
                    case LOG:                   // a*ln(bx), assume b=e/xMax
                        initialParams[0] = yMean;               // don't care, we will do this via regression
                        initialParams[1] = Math.E / (xMax + 1e-100);
                        break;
                    case LOG2:                  // y = a+b*ln(x-c)
                        initialParams[0] = yMean;               // don't care, we will do this via regression
                        initialParams[1] = (yMax - yMin) / (xMax - xMin + 1e-100); // don't care, we will do this via regression
                        initialParams[2] = Math.Min(0.0, -xMin - 0.1 * (xMax - xMin) - 1e-100);
                        break;
                    case RODBARD:               // d+(a-d)/(1+(x/c)^b)
                        initialParams[0] = firsty;  // don't care, we will do this via regression
                        initialParams[1] = 1.0;
                        initialParams[2] = xMin < 0 ? xMin : xMax; //better than xMean;
                        initialParams[3] = lasty;   // don't care, we will do this via regression
                        break;
                    case INV_RODBARD:
                    case RODBARD2: // c*((x-a)/(d-x))^(1/b)
                        initialParams[0] = xMin - 0.1 * (xMax - xMin);
                        initialParams[1] = 1.0;
                        initialParams[2] = yMax;    // don't care, we will do this via regression
                        initialParams[3] = xMax + (xMax - xMin);
                        break;
                    case GAMMA_VARIATE:         // // b*(x-a)^c*exp(-(x-a)/d)
                                                //	First guesses based on following observations (naming the 'x' coordinate 't'):
                                                //	t0 [a] = time of first rise in gamma curve - so use the user specified first x value
                                                //	tmax = t0 + c*d where tmX is the time of the peak of the curve
                                                //	therefore an estimate for c and d is sqrt(tm-t0)
                                                //	K [a] can now be calculated from these estimates
                        initialParams[0] = xMin;
                        double cd = xOfMax - xMin;
                        if (cd < 0.1 * (xMax - xMin)) cd = 0.1 * (xMax - xMin);
                        initialParams[2] = Math.Sqrt(cd);
                        initialParams[3] = Math.Sqrt(cd);
                        initialParams[1] = yMax / (Math.Pow(cd, initialParams[2]) * Math.Exp(-cd / initialParams[3])); // don't care, we will do this via regression
                        break;
                    case GAUSSIAN:                  // a + (b-a)*exp(-(x-c)^2/(2d^2))
                        initialParams[0] = yMin;    //actually don't care, we will do this via regression
                        initialParams[1] = yMax;    //actually don't care, we will do this via regression
                        initialParams[2] = xOfMax;
                        initialParams[3] = 0.39894 * (xMax - xMin) * (yMean - yMin) / (yMax - yMin + 1e-100);
                        break;
                    case GAUSSIAN_NOOFFSET:         // a*exp(-(x-b)^2/(2c^2))
                        initialParams[0] = yMax;    //actually don't care, we will do this via regression
                        initialParams[1] = xOfMax;    //actually don't care, we will do this via regression
                        initialParams[2] = 0.39894 * (xMax - xMin) * yMean / (yMax + 1e-100);
                        break;
                    case CHAPMAN:                   // a*(1-exp(-b*x))^c
                        initialParams[0] = yMax;
                        initialParams[2] = 1.5; // just assuming any reasonable value
                        for (int i = 1; i < numPoints; i++) //find when we reach 50% of maximum
                            if (yData[i] > 0.5 * yMax) {  //approximately (1-exp(-1))^1.5 = 0.5
                                initialParams[1] = 1.0 / xData[i];
                                break;
                            }
                        if (double.IsNaN(initialParams[1]) || initialParams[1] > 1000.0 / xMax) //in case an outlier at the beginning has fooled us
                            initialParams[1] = 10.0 / xMax;
                        break;
                        //no case CUSTOM: here, was done above
                }
            }
            if (!hasInitialParamVariations) {   // estimate initial range for parameters
                for (int i = 0; i < numParams; i++)
                    initialParamVariations[i] = 0.1 * initialParams[i]; //default, should be overridden if it can be zero
                switch (fitType) {
                    case POLY2:
                    case POLY3:
                    case POLY4:
                    case POLY5:
                    case POLY6:
                    case POLY7:
                    case POLY8:
                        double xFactor = 0.5 * Math.Max(Math.Abs(xMax + xMin), (xMax - xMin));
                        initialParamVariations[numParams - 1] = (yMax - yMin) / (Math.Pow(0.5 * (xMax - xMin), numParams - 1) + 1e-100);
                        for (int i = numParams - 2; i >= 0; i--)
                            initialParamVariations[i] = initialParamVariations[i + 1] * xFactor;
                        break;
                    case EXPONENTIAL:        // a*exp(bx)		  a (and c) is calculated by regression
                    case EXP_WITH_OFFSET:    // a*exp(-bx) + c
                    case EXP_RECOVERY:       // a*(1-exp(-bx)) + c
                        initialParamVariations[1] = 0.1 / (xMax - xMin + 1e-100);
                        break;
                    // case CHAPMAN:            // a*(1-exp(-b*x))^c use default (10% of value) for b, c
                    // case POWER:				// ax^b; use default for b
                    // case LOG:				// a*ln(bx); use default for b
                    // case LOG2:				// y = a+b*ln(x-c); use default for c
                    case RODBARD:               // d+(a-d)/(1+(x/c)^b); a and d calculated by regression
                        initialParamVariations[2] = 0.5 * Math.Max((xMax - xMin), Math.Abs(xMean));
                        initialParamVariations[3] = 0.5 * Math.Max(yMax - yMin, Math.Abs(yMax));
                        break;
                    case INV_RODBARD:           // c*((x-a)/(d-x))^(1/b); c calculated by regression
                        initialParamVariations[0] = 0.01 * Math.Max(xMax - xMin, Math.Abs(xMax));
                        initialParamVariations[2] = 0.1 * Math.Max(yMax - yMin, Math.Abs(yMax));
                        initialParamVariations[3] = 0.1 * Math.Max((xMax - xMin), Math.Abs(xMean));
                        break;
                    case GAMMA_VARIATE:         // // b*(x-a)^c*exp(-(x-a)/d); b calculated by regression
                                                //	First guesses based on following observations:
                                                //	t0 [b] = time of first rise in gamma curve - so use the user specified first limit
                                                //	tm = t0 + a*B [c*d] where tm is the time of the peak of the curve
                                                //	therefore an estimate for a and B is sqrt(tm-t0)
                                                //	K [a] can now be calculated from these estimates
                        initialParamVariations[0] = 0.1 * Math.Max(yMax - yMin, Math.Abs(yMax));
                        double ab = xOfMax - firstx + 0.1 * (xMax - xMin);
                        initialParamVariations[2] = 0.1 * Math.Sqrt(ab);
                        initialParamVariations[3] = 0.1 * Math.Sqrt(ab);
                        break;
                    case GAUSSIAN:              // a + (b-a)*exp(-(x-c)^2/(2d^2)); a,b calculated by regression
                        initialParamVariations[2] = 0.2 * initialParams[3]; //(and default for d)
                        break;
                    case GAUSSIAN_NOOFFSET:     // a*exp(-(x-b)^2/(2c^2))
                        initialParamVariations[1] = 0.2 * initialParams[2]; //(and default for c)
                        break;
                }
            }
            return true;
        }

        /** Set multiplyParams and offsetParam for built-in functions. This allows us to use linear
         *	regression, reducing the number of parameters used by the io.binroot.regression.Minimizer by up to 2, and
         *	improving the speed and success rate of the minimization process */
        private void getOffsetAndFactorParams() {
            offsetParam = -1;
            factorParam = -1;
            hasSlopeParam = false;
            switch (fitType) {
                case STRAIGHT_LINE:
                case POLY2:
                case POLY3:
                case POLY4:
                case POLY5:
                case POLY6:
                case POLY7:
                case POLY8:
                    offsetParam = 0;
                    factorParam = 1;
                    hasSlopeParam = true;
                    break;
                case EXPONENTIAL:           // a*exp(bx)
                    factorParam = 0;
                    break;
                case EXP_WITH_OFFSET:       // a*exp(-bx) + c
                case EXP_RECOVERY:          // a*(1-exp(-bx)) + c
                    offsetParam = 2;
                    factorParam = 0;
                    break;
                case EXP_RECOVERY_NOOFFSET: // a*(1-exp(-bx))
                    factorParam = 0;
                    break;
                case CHAPMAN:               // a*(1-exp(-b*x))^c
                    factorParam = 0;
                    break;
                case POWER:                 // ax^b
                    factorParam = 0;
                    break;
                case LOG:                   // a*ln(bx)
                    factorParam = 0;
                    break;
                case LOG2:                  // y = a+b*ln(x-c)
                    offsetParam = 0;
                    factorParam = 1;
                    break;
                case RODBARD_INTERNAL:      // d+a/(1+(x/c)^b)
                    offsetParam = 3;
                    factorParam = 0;
                    break;
                case INV_RODBARD:           // c*((x-a)/(d-x))^(1/b)
                    factorParam = 2;
                    break;
                case GAMMA_VARIATE:         // b*(x-a)^c*exp(-(x-a)/d)
                    factorParam = 1;
                    break;
                case GAUSSIAN_INTERNAL:     // a + b*exp(-(x-c)^2/(2d^2))
                    offsetParam = 0;
                    factorParam = 1;
                    break;
                case GAUSSIAN_NOOFFSET:     // a*exp(-(x-b)^2/(2c^2))
                    factorParam = 0;
                    break;
            }
            numRegressionParams = 0;
            if (offsetParam >= 0) numRegressionParams++;
            if (factorParam >= 0) numRegressionParams++;
        }


        /** calculates the sum of y and y^2 */
        private void calculateSumYandY2() {
            sumY = 0.0; sumY2 = 0.0;
            for (int i = 0; i < numPoints; i++) {
                double y = yData[i];
                sumY += y;
                sumY2 += y * y;
            }
        }

        /** returns whether this a fit type that acutally fits modified data with a modified function */
        private bool isModifiedFitType(int fitType) {
            return fitType == POWER_REGRESSION || fitType == EXP_REGRESSION || fitType == RODBARD2 ||
                    fitType == GAUSSIAN;
        }

        /** For fits don't use the original data, prepare modified data and fit type.
         *	Returns false if the data are incompatible with the fit type
         *	In that case, 'errorString' is set to a message explaining the problem
         */
        private bool prepareModifiedFitType(int fitType) {
            if (fitType == GAUSSIAN) {
                this.fitType = GAUSSIAN_INTERNAL;   // different definition of parameters for using regression
                return true;
            } else if (fitType == RODBARD) {
                this.fitType = RODBARD_INTERNAL;    // different definition of parameters for using regression
                return true;
            } else if (fitType == POWER_REGRESSION || fitType == EXP_REGRESSION) {
                if (fitType == POWER_REGRESSION) {
                    xDataSave = xData;
                    xData = new double[numPoints];
                }
                yDataSave = yData;
                yData = new double[numPoints];
                ySign = 0;
                numPoints = 0;  // we may have lower number of points if there is a (0,0) point that we ignore
                for (int i = 0; i < xData.Length; i++) {
                    double y = yDataSave[i];
                    if (fitType == POWER_REGRESSION) {
                        double x = xDataSave[i];
                        if (x == 0 && y == 0) {
                            restrictPower = true;
                            continue;     // ignore (0,0) point in power law
                        }
                        if (x <= 0) {
                            errorString = "Cannot fit x<=0";
                            return false;
                        }
                        xData[numPoints] = Math.Log(x);
                    }
                    if (ySign == 0) ySign = Math.Sign(y); //if unknown, determine whether y data are positive or negative
                    if (y * ySign <= 0) {
                        errorString = "Cannot fit y=0 or mixture of y>0, y<0";
                        return false;
                    }
                    yData[numPoints] = Math.Log(y * ySign);
                    numPoints++;
                }
                this.fitType = STRAIGHT_LINE;
            } else if (fitType == RODBARD2) { // 'Rodbard (NIH Image)' fit is inverse Rodbard function, by fitting x(y); for compatibility with NIH Image
                xDataSave = xData;
                yDataSave = yData;
                xData = yDataSave;  //swap
                yData = xDataSave;
                this.fitType = RODBARD_INTERNAL;
            }
            return true;
        }

        /** Get correct params1 and revert xData, yData if fit has been done via another function */
        private void postProcessModifiedFitType(int fitType) {
            if (fitType == POWER_REGRESSION || fitType == EXP_REGRESSION)   // ln y = ln (a*x^b) = ln a + b ln x
                finalParams[0] = ySign * Math.Exp(finalParams[0]);          //or: ln (+,-)y = ln ((+,-)a*exp(bx)) = ln (+,-)a + bx
            if (fitType == GAUSSIAN)                            // a + b exp(-...) to  a + (b-a)*exp(-...)
                finalParams[1] += finalParams[0];
            else if (fitType == RODBARD || fitType == RODBARD2) //d+a/(1+(x/c)^b) to d+(a-d)/(1+(x/c)^b)
                finalParams[0] += finalParams[3];

            if (xDataSave != null) {
                xData = xDataSave;
                numPoints = xData.Length;
            }
            if (yDataSave != null) yData = yDataSave;
            this.fitType = fitType;
        }

        private double sqr(double d) { return d * d; }

        /**
         * Gets index of highest value in an array.
         *
         * @param			   array the array.
         * @return			   Index of highest value.
         */
        public static int getMax(double[] array) {
            double max = array[0];
            int index = 0;
            for (int i = 1; i < array.Length; i++) {
                if (max < array[i]) {
                    max = array[i];
                    index = i;
                }
            }
            return index;
        }

    }
}
