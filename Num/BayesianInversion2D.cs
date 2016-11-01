/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;

namespace MaxQuant.Num {
	public class BayesianInversion2D {
		private readonly double[] xout;
		private readonly double[] yout;
		private readonly double[,] zout;
		private readonly double[,] forward;
		private readonly double[,] reverse;

		public BayesianInversion2D(double[] xdata, double[] ydata, bool[] correct) : this(xdata, ydata, correct, false) {
		}

		public BayesianInversion2D(double[] xdata, double[] ydata, bool[] correct, bool debug) {
			Invert(xdata, ydata, correct, out xout, out yout, out zout, debug, out forward, out reverse);
		}

		public double GetValue(double x, double y) {
			return GetValue(x, y, xout, yout, zout);
		}

		public double GetForwardHist(double x, double y) {
			return GetValue(x, y, xout, yout, forward);
		}

		public double GetReverseHist(double x, double y) {
			return GetValue(x, y, xout, yout, reverse);
		}

		private static void Invert(double[] xdata, double[] ydata, bool[] correct, out double[] xout, out double[] yout,
		                           out double[,] zout, bool debug, out double[,] forwardOut, out double[,] reverseOut) {
			int n = correct.Length;
			int ntrue = 0;
			for (int i = 0; i < n; i++) {
				if (correct[i]) {
					ntrue++;
				}
			}
			int nfalse = n - ntrue;
			if (ntrue < 3 || nfalse < 3) {
				xout = null;
				yout = null;
				zout = null;
				forwardOut = null;
				reverseOut = null;
				return;
			}
			double falseP = nfalse / (double) n;
			double[,] falseData = new double[nfalse,2];
			double[,] trueData = new double[ntrue,2];
			int iFalse = 0;
			int iTrue = 0;
			for (int i = 0; i < n; i++) {
				if (correct[i]) {
					trueData[iTrue, 0] = xdata[i];
					trueData[iTrue, 1] = ydata[i];
					iTrue++;
				} else {
					falseData[iFalse, 0] = xdata[i];
					falseData[iFalse, 1] = ydata[i];
					iFalse++;
				}
			}
			Invert(falseData, trueData, falseP, out xout, out yout, out zout, out forwardOut, out reverseOut, debug);
		}

		private static double GetValue(double x, double y, double[] xvals, double[] yvals, double[,] zvals) {
			double z = GetValueImpl(x, y, xvals, yvals, zvals);
			if(double.IsNaN(z)){
				return double.NaN;
			}
			//return Math.Min(1, z);
			return z;
		}

		private static double GetValueImpl(double x, double y, double[] xvals, double[] yvals, double[,] zvals) {
			if (xvals == null || yvals == null) {
				return double.NaN;
			}
			if (xvals.Length == 0 || yvals.Length == 0) {
				return double.NaN;
			}
			if (x <= xvals[0]) {
				return InterpolateExactX(0, y, yvals, zvals);
			}
			if (x >= xvals[xvals.Length - 1]) {
				return InterpolateExactX(xvals.Length - 1, y, yvals, zvals);
			}
			int ax = Array.BinarySearch(xvals, x);
			if (ax >= 0) {
				return InterpolateExactX(ax, y, yvals, zvals);
			}
			if (y <= yvals[0]) {
				return InterpolateExactY(0, x, xvals, zvals);
			}
			if (y >= yvals[yvals.Length - 1]) {
				return InterpolateExactY(yvals.Length - 1, x, xvals, zvals);
			}
			int ay = Array.BinarySearch(yvals, y);
			if (ay >= 0) {
				return InterpolateExactY(ay, x, xvals, zvals);
			}
			int i1 = -2 - ax;
			int i2 = i1 + 1;
			int j1 = -2 - ay;
			int j2 = j1 + 1;
			if (i1 < 0 || j1 < 0) {
				return double.NaN;
			}
			if (i2 >= xvals.Length || j2 >= yvals.Length) {
				return double.NaN;
			}
			double x1 = xvals[i1];
			double x2 = xvals[i2];
			double y1 = yvals[j1];
			double y2 = yvals[j2];
			double z11 = zvals[i1, j1];
			double z12 = zvals[i1, j2];
			double z21 = zvals[i2, j1];
			double z22 = zvals[i2, j2];
			double w1 = (z11 * (x2 - x) + z21 * (x - x1)) / (x2 - x1);
			double w2 = (z12 * (x2 - x) + z22 * (x - x1)) / (x2 - x1);
			return (w1 * (y2 - y) + w2 * (y - y1)) / (y2 - y1);
		}

		private static void Invert(double[,] falseData, double[,] trueData, double falseP, out double[] xRes,
		                           out double[] yRes, out double[,] zRes, out double[,] forwardOut, out double[,] reverseOut,
		                           bool debug) {
			double xmin = double.MaxValue;
			double xmax = -double.MaxValue;
			double ymin = double.MaxValue;
			double ymax = -double.MaxValue;
			for (int i = 0; i < falseData.GetLength(0); i++) {
				double xx = falseData[i, 0];
				double yy = falseData[i, 1];
				if (xx < xmin) {
					xmin = xx;
				}
				if (xx > xmax) {
					xmax = xx;
				}
				if (yy < ymin) {
					ymin = yy;
				}
				if (yy > ymax) {
					ymax = yy;
				}
			}
			for (int i = 0; i < trueData.GetLength(0); i++) {
				double xx = trueData[i, 0];
				double yy = trueData[i, 1];
				if (xx < xmin) {
					xmin = xx;
				}
				if (xx > xmax) {
					xmax = xx;
				}
				if (yy < ymin) {
					ymin = yy;
				}
				if (yy > ymax) {
					ymax = yy;
				}
			}
			double dx = xmax - xmin;
			xmin -= 0.1 * dx;
			xmax += 0.1 * dx;
			double dy = ymax - ymin;
			ymin -= 0.1 * dy;
			ymax += 0.1 * dy;
			double[] falseX;
			double[] falseY;
			double[,] falseZ;
			int n = trueData.GetLength(0);
			double[,] cov = NumUtil.CalcCovariance(trueData);
			double fact = Math.Pow(n, 1.0 / 6.0);
			double[,] hinv = null;
			try {
				hinv = NumUtil.ApplyFunction(cov, delegate(double w) { return fact / Math.Sqrt(w); });
			} catch {
			}
			if (hinv == null || !IsValidMatrix(hinv)) {
				xRes = null;
				yRes = null;
				zRes = null;
				forwardOut = null;
				reverseOut = null;
				return;
			}
			EstimateBivariateDensity(falseData, out falseX, out falseY, out falseZ, xmin, xmax, ymin, ymax, hinv);
			double[] trueX;
			double[] trueY;
			double[,] trueZ;
			EstimateBivariateDensity(trueData, out trueX, out trueY, out trueZ, xmin, xmax, ymin, ymax, hinv);
			double[] x = UnifySupport(falseX, trueX);
			double[] y = UnifySupport(falseY, trueY);
			falseZ = Interpolate(x, y, falseX, falseY, falseZ);
			trueZ = Interpolate(x, y, trueX, trueY, trueZ);
			double[,] inverse = new double[x.Length,y.Length];
			for (int i = 0; i < x.Length; i++) {
				for (int j = 0; j < y.Length; j++) {
					inverse[i, j] = falseZ[i, j] <= 0 ? double.Epsilon : Math.Max((falseZ[i, j] * 0.5) / trueZ[i, j], double.Epsilon);
				}
			}
			for (int j = 0; j < y.Length; j++) {
				double maxVal = double.MinValue;
				int maxInd = -1;
				for (int i = 0; i < x.Length; i++) {
					if(inverse[i,j] > maxVal) {
						maxVal = inverse[i, j];
						maxInd = i;
					}
				}
				for (int i = 0; i < maxInd; i++) {
					inverse[i, j] = maxVal;
				}
			}
			xRes = x;
			yRes = y;
			zRes = inverse;
			if (debug) {
				forwardOut = trueZ;
				reverseOut = falseZ;
			} else {
				forwardOut = null;
				reverseOut = null;
			}
		}

		private static bool IsValidMatrix(double[,] x) {
			for (int i = 0; i < x.GetLength(0); i++) {
				for (int j = 0; j < x.GetLength(1); j++) {
					double y = x[i, j];
					if (double.IsNaN(y) || double.IsInfinity(y)) {
						return false;
					}
				}
			}
			return true;
		}

		private static double[] UnifySupport(double[] x1, double[] x2) {
			if (x1.Length == 0) {
				return x2;
			}
			if (x2.Length == 0) {
				return x1;
			}
			double avDiff1 = 0;
			for (int i = 0; i < x1.Length - 1; i++) {
				avDiff1 += x1[i + 1] - x1[i];
			}
			avDiff1 /= x1.Length - 1;
			double avDiff2 = 0;
			for (int i = 0; i < x2.Length - 1; i++) {
				avDiff2 += x2[i + 1] - x2[i];
			}
			avDiff2 /= x2.Length - 1;
			double diff = Math.Min(avDiff1, avDiff2);
			double min = Math.Min(x1[0], x2[0]);
			double max = Math.Max(x1[x1.Length - 1], x2[x2.Length - 1]);
			List<double> xvals = new List<double>();
			for (double val = min; val <= max; val += diff) {
				xvals.Add(val);
			}
			return xvals.ToArray();
		}

		private static double[,] Interpolate(double[] newX, double[] newY, double[] xvals, double[] yvals, double[,] zvals) {
			double[,] newZ = new double[newX.Length,newY.Length];
			for (int i = 0; i < newX.Length; i++) {
				for (int j = 0; j < newY.Length; j++) {
					newZ[i, j] = GetValue(newX[i], newY[j], xvals, yvals, zvals);
				}
			}
			return newZ;
		}

		private static double InterpolateExactY(int yind, double x, double[] xvals, double[,] zvals) {
			if (x <= xvals[0]) {
				return zvals[0, yind];
			}
			if (x >= xvals[xvals.Length - 1]) {
				return zvals[xvals.Length - 1, yind];
			}
			int ax = Array.BinarySearch(xvals, x);
			if (ax >= 0) {
				return zvals[ax, yind];
			}
			int i1 = -2 - ax;
			int i2 = i1 + 1;
			double x1 = xvals[i1];
			double x2 = xvals[i2];
			double w1 = zvals[i1, yind];
			double w2 = zvals[i2, yind];
			return (w1 * (x2 - x) + w2 * (x - x1)) / (x2 - x1);
		}

		private static double InterpolateExactX(int xind, double y, double[] yvals, double[,] zvals) {
			if (y <= yvals[0]) {
				return zvals[xind, 0];
			}
			if (y >= yvals[yvals.Length - 1]) {
				return zvals[xind, yvals.Length - 1];
			}
			int ay = Array.BinarySearch(yvals, y);
			if (ay >= 0) {
				return zvals[xind, ay];
			}
			int j1 = -2 - ay;
			int j2 = j1 + 1;
			double y1 = yvals[j1];
			double y2 = yvals[j2];
			double w1 = zvals[xind, j1];
			double w2 = zvals[xind, j2];
			return (w1 * (y2 - y) + w2 * (y - y1)) / (y2 - y1);
		}

		private static void EstimateBivariateDensity(double[,] data, out double[] xvals, out double[] yvals,
		                                             out double[,] zvals, double xmin, double xmax, double ymin, double ymax,
		                                             double[,] hinv) {
			double bandWidthX = 1.0 / hinv[0, 0];
			double bandWidthY = 1.0 / hinv[1, 1];
			List<double> xv = new List<double>();
			for (double val = xmin; val <= xmax; val += bandWidthX) {
				xv.Add(val);
			}
			xvals = xv.ToArray();
			List<double> yv = new List<double>();
			for (double val = ymin; val <= ymax; val += bandWidthY) {
				yv.Add(val);
			}
			yvals = yv.ToArray();
			zvals = new double[xvals.Length,yvals.Length];
			for (int i = 0; i < xvals.Length; i++) {
				for (int j = 0; j < yvals.Length; j++) {
					zvals[i, j] = EstimateDensity(xvals[i], yvals[j], data, hinv);
				}
			}
		}

		private static double EstimateDensity(double x, double y, double[,] data, double[,] hinv) {
			double result = 0;
			for (int i = 0; i < data.GetLength(0); i++) {
				double[] w = new double[] {x - data[i, 0], y - data[i, 1]};
				double[] b = NumUtil.MatrixTimesVector(hinv, w);
				result += NumUtil.StandardGaussian(b);
			}
			result *= NumUtil.Determinant2x2(hinv) / data.Length;
			return result;
		}

		public bool IsValid() {
			return zout != null;
		}
	}
}