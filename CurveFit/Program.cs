using System;

namespace CurveFit {
    class Program {
        static void Main(string[] args) {
            double[] xs = { 1, 2, 3, 4, 5, 6, 7, 14, 30 };
            double[] ys = { 0.228, 0.1638, 0.1402, 0.1306, 0.1229, 0.1149, 0.1104, 0.0954, 0.0687 };
            var curveFitter = new CurveFitter(xs, ys);
            curveFitter.doFit(CurveFitter.POWER);
            var params1 = curveFitter.getParams();

            var a = params1[0];
            var b = params1[1];
            Console.WriteLine("y = {0} * pow(x, {1})", a, b);
            for (int i = 0; i < xs.Length; i++) {
                var x = xs[i];
                var ye = ys[i];
                var ya = a * Math.Pow(x, b);
                Console.WriteLine("{0}, {1}, {2:f4}, {3:f4}", x, ye, ya, ya - ye);
            }
            foreach(var x in new int[] { 60, 90, 120, 365}) {
                var y = a * Math.Pow(x, b);
                Console.WriteLine("{0}, {1:f4}", x, y);
            }
            Console.ReadLine();
        }
    }
}
