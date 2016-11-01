/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;

namespace MaxQuant.Num.Test {
	public class OneSampleTTest : OneSampleTest {
		public override void Test(double[] data, double mean, out double statistic, out double bothTails,
		                          out double leftTail, out double rightTail) {
			TestImpl(data, mean, out statistic, out bothTails, out leftTail, out rightTail);
		}

		private static void TestImpl(double[] x, double mean, out double stat, out double bothtails,
		                             out double lefttail, out double righttail) {
			int n = x.Length;
			int i;
			stat = 0;

			if (n <= 1) {
				bothtails = 1.0;
				lefttail = 1.0;
				righttail = 1.0;
				return;
			}

			//
			// Mean
			//
			double xmean = 0;
			for (i = 0; i <= n - 1; i++) {
				xmean = xmean + x[i];
			}
			xmean = xmean / n;

			//
			// Variance (using corrected two-pass algorithm)
			//
			double xstddev = 0;
			if (n != 1) {
				double v1 = 0;
				for (i = 0; i <= n - 1; i++) {
					v1 = v1 + (x[i] - xmean) * (x[i] - xmean);
				}
				double v2 = 0;
				for (i = 0; i <= n - 1; i++) {
					v2 = v2 + (x[i] - xmean);
				}
				v2 = v2 * v2 / n;
				double xvariance = (v1 - v2) / (n - 1);
				if (xvariance < 0) {
					xvariance = 0;
				}
				xstddev = Math.Sqrt(xvariance);
			}
			if (xstddev == 0) {
				bothtails = 1.0;
				lefttail = 1.0;
				righttail = 1.0;
				return;
			}
			stat = (xmean - mean) / (xstddev / Math.Sqrt(n));
			double df = n - 1;
			try {
				double p;
				if (stat > 0) {
					p = 1 - 0.5 * NumericalRecipes.Betai(df / 2, 0.5, df / (df + stat * stat));
				} else {
					p = 0.5 * NumericalRecipes.Betai(df / 2, 0.5, df / (df + stat * stat));
				}
				bothtails = 2 * Math.Min(p, 1 - p);
				lefttail = p;
				righttail = 1 - p;
			} catch (Exception) {
				bothtails = 1.0;
				lefttail = 1.0;
				righttail = 1.0;
				return;
			}
		}
	}
}