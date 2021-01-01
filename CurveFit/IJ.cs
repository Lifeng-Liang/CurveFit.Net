using System;

namespace CurveFit {
    public class IJ {
        /** Converts a number to a rounded formatted string.
         The 'decimalPlaces' argument specifies the number of
         digits to the right of the decimal point (0-9). Uses
         scientific notation if 'decimalPlaces is negative. */
        public static string d2s(double n, int decimalPlaces) {
            return string.Format("{0:f" + decimalPlaces + "}", n);
        }

        /** Converts a number to a rounded formatted string.
         * The 'significantDigits' argument specifies the minimum number
         * of significant digits, which is also the preferred number of
         * digits behind the decimal. Fewer decimals are shown if the
         * number would have more than 'maxDigits'.
         * Exponential notation is used if more than 'maxDigits' would be needed.
         */
        public static string d2s(double x, int significantDigits, int maxDigits) {
            double log10 = Math.Log10(Math.Abs(x));
            double roundErrorAtMax = 0.223 * Math.Pow(10, -maxDigits);
            int magnitude = (int)Math.Ceiling(log10 + roundErrorAtMax);
            int decimals = x == 0 ? 0 : maxDigits - magnitude;
            if (decimals < 0 || magnitude < significantDigits + 1 - maxDigits)
                return IJ.d2s(x, -significantDigits); // exp notation for large and small numbers
            else {
                if (decimals > significantDigits)
                    decimals = Math.Max(significantDigits, decimals - maxDigits + significantDigits);
                return IJ.d2s(x, decimals);
            }
        }
    }
}
