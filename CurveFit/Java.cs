using System;
using System.Collections.Generic;

namespace CurveFit {
    public static class Java {
        public static List<double[]> NewDoubleArray2(int x, int y) {
            var list = new List<double[]>();
            for (int i = 0; i < x; i++) {
                list.Add(new double[y]);
            }
            return list;
        }

        public static long SystemCurrentTimeMillis() {
            var currentTicks = DateTime.Now.Ticks;
            var dtFrom = new DateTime(1970, 1, 1, 0, 0, 0, 0);
            var currentMillis = (currentTicks - dtFrom.Ticks) / 10000;
            return currentMillis;
        }

        public static void DoubleArrayFill(double[] array, double v) {
            for (int i = 0; i < array.Length; i++) {
                array[i] = v;
            }
        }
    }
}
