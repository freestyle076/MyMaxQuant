/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using MaxQuant.Util;

namespace MaxQuant.Num {
	public delegate double Func(double x);
	public delegate double Funcs2(double x, double[] a);

	public static class NumUtil {

		public static double log10 = Math.Log(10);

		private const double EPS = 1.0e-7;

		private const double mass0 = 445.120025;

		public static double RelativeCorrection(double u, double[,] mzCalibrationPar, double[,] intensityCalibrationPar,
		                                        double intens, int massRange) {
			if (mzCalibrationPar == null) {
				return 0;
			}
			u = 100 / Math.Sqrt(u) - 100 / Math.Sqrt(mass0);
			double w = Math.Log(intens) - Math.Log(1e7);
			double result = 0;
			double fact = u;
			for (int i = 0; i < mzCalibrationPar.GetLength(1); i++) {
				result += fact * mzCalibrationPar[massRange, i];
				fact *= u;
			}
			fact = w;
			for (int i = 0; i < intensityCalibrationPar.GetLength(1); i++) {
				result += fact * intensityCalibrationPar[massRange, i];
				fact *= w;
			}
			return result;
		}

		public static void LinFit2(double[] x, double[] y, double[] a, LfitFunc funcs) {
			double chisq;
			LinFit2(x, y, null, a, out chisq, funcs);
		}

		public static void LinFit2(double[] x, double[] y, double[] sig, double[] a,
		                           out double chisq, LfitFunc funcs) {
			double[,] covar;
			if (sig == null) {
				sig = new double[x.Length];
				for (int i = 0; i < sig.Length; i++) {
					sig[i] = 1E-2;
				}
			}
			NumericalRecipes.Lfit(x, y, sig, a, out covar, out chisq, funcs);
		}

		public static int[][] GetPartitions(int n, int ngroups) {
			List<int[]> partitions = new List<int[]>();
			Partition(new TmpPartition(n), partitions, ngroups);
			return partitions.ToArray();
		}

		public static double Multinomial(int n, int[] partition) {
			return Math.Exp(LnMultinomial(n, partition));
		}

		public static double LnMultinomial(int a, int[] bs) {
			double result = NumericalRecipes.Gammln(a + 1);
			for (int i = 0; i < bs.Length; i++) {
				result -= NumericalRecipes.Gammln(bs[i] + 1);
			}
			return result;
		}

		private static void Partition(TmpPartition x, ICollection<int[]> allPartitions, int len) {
			if (x.remainder == 0 && x.partition.Count == len) {
				allPartitions.Add(x.partition.ToArray());
				return;
			}
			if (x.partition.Count == len) {
				return;
			}
			for (int i = 0; i <= x.remainder; i++) {
				Partition(x.add(i), allPartitions, len);
			}
		}

		public static void FitOrigin(double[] x, double[] y, double[] sig,
		                             out double b, out double sigb, out double chi2, out double q) {
			int ndata = x.Length;
			if (x.Length != y.Length) {
				throw new Exception("x.Length != y.Length: " + x.Length + " " + y.Length);
			}
			if (sig != null && x.Length != sig.Length) {
				throw new Exception("x.Length != sig.Length: " + x.Length + " " + sig.Length);
			}
			double t;
			double st2 = 0.0;
			b = 0.0;
			if (sig != null) {
				for (int i = 0; i < ndata; i++) {
					t = (x[i]) / sig[i];
					st2 += t * t;
					b += t * y[i] / sig[i];
				}
			} else {
				for (int i = 0; i < ndata; i++) {
					t = x[i];
					st2 += t * t;
					b += t * y[i];
				}
			}
			b /= st2;
			sigb = Math.Sqrt(1.0 / st2);
			chi2 = 0.0;
			if (sig == null) {
				for (int i = 0; i < ndata; i++) {
					double tmp = y[i] - b * x[i];
					chi2 += tmp * tmp;
				}
				q = 1.0;
				double sigdat = Math.Sqrt((chi2) / (ndata - 2));
				sigb *= sigdat;
			} else {
				for (int i = 0; i < ndata; i++) {
					double tmp = (y[i] - b * x[i]) / sig[i];
					chi2 += tmp * tmp;
				}
				q = NumericalRecipes.Gammq(0.5 * (ndata - 1), 0.5 * (chi2));
			}
		}

		public static double MedfitOrigin(double[] x, double[] y, out double abdev) {
			int ndata = x.Length;
			int j;
			double sx = 0.0, sy = 0.0, sxy = 0.0, sxx = 0.0, chisq = 0.0;
			int ndatat = ndata;
			double[] xt = x;
			double[] yt = y;
			for (j = 0; j < ndata; j++) {
				sx += x[j];
				sy += y[j];
				sxy += x[j] * y[j];
				sxx += x[j] * x[j];
			}
			double del = ndata * sxx - sx * sx;
			double bb = (ndata * sxy - sx * sy) / del;
			for (j = 0; j < ndata; j++) {
				double temp = y[j] - (bb * x[j]);
				chisq += (temp * temp);
			}
			double sigb = Math.Sqrt(chisq / del);
			double b1 = bb;
			double abdevt = 0;
			double f1 = Rofunc(b1, ndatat, xt, yt, ref abdevt);
			double sign = f1 >= 0.0 ? Math.Abs(3.0 * sigb) : -Math.Abs(3.0 * sigb);
			double b2 = bb + sign;
			double f2 = Rofunc(b2, ndatat, xt, yt, ref abdevt);
			while (f1 * f2 > 0.0) {
				bb = 2.0 * b2 - b1;
				b1 = b2;
				f1 = f2;
				b2 = bb;
				f2 = Rofunc(b2, ndatat, xt, yt, ref abdevt);
			}
			sigb = 0.01 * sigb;
			while (Math.Abs(b2 - b1) > sigb) {
				bb = 0.5 * (b1 + b2);
				if (bb == b1 || bb == b2) {
					break;
				}
				double f = Rofunc(bb, ndatat, xt, yt, ref abdevt);
				if (f * f1 >= 0.0) {
					f1 = f;
					b1 = bb;
				} else {
					b2 = bb;
				}
			}
			abdev = abdevt / ndata;
			return bb;
		}

		public static double Rofunc(double b, int ndatat, double[] xt, double[] yt, ref double abdevt) {
			int j;
			double sum = 0.0;
			double[] arr = new double[ndatat];
			for (j = 0; j < ndatat; j++) {
				arr[j] = yt[j] - b * xt[j];
			}
			abdevt = 0.0;
			for (j = 0; j < ndatat; j++) {
				double d = yt[j] - (b * xt[j]);
				abdevt += Math.Abs(d);
				if (yt[j] != 0.0) {
					d /= Math.Abs(yt[j]);
				}
				if (Math.Abs(d) > EPS) {
					sum += (d >= 0.0 ? xt[j] : -xt[j]);
				}
			}
			return sum;
		}

		private static double Dyda(double x, double[] a, Funcs2 func, int ind, double epsilon) {
			double[] a1 = new double[a.Length];
			double[] a2 = new double[a.Length];
			for (int i = 0; i < a.Length; i++) {
				a1[i] = a[i];
				a2[i] = a[i];
			}
			a1[ind] += epsilon / 2.0;
			a2[ind] -= epsilon / 2.0;
			return (func(x, a1) - func(x, a2)) / epsilon;
		}

		public static void FitNonlin(double[] x, double[] y, double[] sig, double[] a,
		                             out double chisq, Funcs2 func) {
			MrqminFunc f = delegate(double x1, double[] a1, out double y1, double[] dyda, int na) {
				y1 = func(x1, a1);
				for (int i = 0; i < na; i++) {
					dyda[i] = Dyda(x1, a1, func, i, 1e-6);
				}
			};
			FitNonlin(x, y, sig, a, out chisq, f);
		}

		public static void FitNonlin(double[] x, double[] y, double[] sig, double[] a,
		                             out double chisq, MrqminFunc func) {
			int ndata = x.Length;
			if (sig == null) {
				sig = new double[ndata];
				for (int i = 0; i < sig.Length; i++) {
					sig[i] = 1;
				}
			}
			int ma = a.Length;
			double[,] covar = new double[ma,ma];
			double[,] alpha = new double[ma,ma];
			int[] ia = new int[ma];
			for (int i = 0; i < ma; i++) {
				ia[i] = 1;
			}
			double alamda = -1;
			double ochisq = 0;
			double[,] oneda = null;
			int mfit = 0;
			double[] atry = null;
			double[] beta = null;
			double[] da = null;
			NumericalRecipes.Mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, out chisq, func, ref alamda, ref ochisq,
			                        ref oneda, ref mfit, ref atry, ref beta, ref da);
			int count1 = 0;
			while (alamda > 1e-20 && alamda < 1e20 && count1 < 100) {
				NumericalRecipes.Mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, out chisq, func, ref alamda, ref ochisq,
				                        ref oneda, ref mfit, ref atry, ref beta, ref da);
				count1++;
			}
			alamda = 0;
			NumericalRecipes.Mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, out chisq, func, ref alamda, ref ochisq,
			                        ref oneda, ref mfit, ref atry, ref beta, ref da);
		}

		public static int[][] GetCombinations(int n, int k, int max, out bool incomplete) {
			List<int[]> result = new List<int[]>();
			Combination comb = new Combination(n, k);
			result.Add(comb.Data);
			incomplete = false;
			int count1 = 1;
			while ((comb = comb.Successor) != null) {
				result.Add(comb.Data);
				count1++;
				if (count1 >= max) {
					incomplete = true;
					break;
				}
			}
			return result.ToArray();
		}

		public static double[,] CalcCovariance(double[,] data) {
			int n = data.GetLength(0);
			int p = data.GetLength(1);
			double[] means = new double[p];
			for (int i = 0; i < p; i++) {
				for (int j = 0; j < n; j++) {
					means[i] += data[j, i];
				}
				means[i] /= n;
			}
			double[,] cov = new double[p, p];
			for (int i = 0; i < p; i++) {
				for (int j = 0; j <= i; j++) {
					for (int k = 0; k < n; k++) {
						cov[i, j] += (data[k, i] - means[i]) * (data[k, j] - means[j]);
					}
					cov[i, j] /= n;
					cov[j, i] = cov[i, j];
				}
			}
			return cov;
		}

		public static double[,] ApplyFunction(double[,] m, Func func) {
			int n = m.GetLength(0);
			double[,] v;
			double[] e = DiagonalizeSymmMatrix(m, out v);
			for (int i = 0; i < n; i++) {
				e[i] = func(e[i]);
			}
			double[,] result = new double[n, n];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					for (int k = 0; k < n; k++) {
						result[i, j] += v[i, k] * e[k] * v[j, k];
					}
				}
			}
			return result;
		}

		public static double[] DiagonalizeSymmMatrix(double[,] m, out double[,] evec) {
			double[] d;
			double[] e;
			evec = (double[,])m.Clone();
			NumericalRecipes.Tred2(evec, out d, out e);
			NumericalRecipes.Tqli(d, e, evec);
			return d;
		}

		public static double[] MatrixTimesVector(double[,] x, double[] v) {
			double[] result = new double[x.GetLength(0)];
			for (int i = 0; i < x.GetLength(0); i++) {
				for (int k = 0; k < x.GetLength(1); k++) {
					result[i] += x[i, k] * v[k];
				}
			}
			return result;
		}

		public static double StandardGaussian(double[] x) {
			double sum = 0;
			for (int i = 0; i < x.Length; i++) {
				sum += x[i] * x[i];
			}
			return Math.Exp(-0.5 * sum) / Math.Pow(2 * Math.PI, 0.5 * x.Length);
		}

		public static double Determinant2x2(double[,] m) {
			return m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0];
		}

		public static double[] MovingBoxPlot(double[] values, double[] controlValues,
									 int nbins) {
			return MovingBoxPlot(values, controlValues,
								 nbins, SignificanceType.both);
		}

		public static double[] MovingBoxPlot(double[] values, double[] controlValues,
											 int nbins, SignificanceType type) {
			double[] lowerQuart;
			double[] median;
			double[] upperQuart;
			double[] binBoundaries;
			return MovingBoxPlot(values, controlValues, nbins, out lowerQuart, out median,
								 out upperQuart, out binBoundaries, type);
		}

		public static readonly int minBinsize = 500;

		public static double[] MovingBoxPlot(double[] values, double[] controlValues, int nbins, out double[] lowerQuart,
									 out double[] median, out double[] upperQuart, out double[] binBoundaries,
									 SignificanceType type) {
			int n = values.Length;
			if (n == 0) {
				lowerQuart = new double[0];
				median = new double[0];
				upperQuart = new double[0];
				binBoundaries = new double[0];
				return new double[0];
			}
			nbins = nbins < 0 ? Math.Max(1, (int)Math.Round(n / (double)minBinsize)) : Math.Min(nbins, 1 + n / 4);
			int[] pos = new int[nbins + 1];
			for (int i = 0; i < nbins + 1; i++) {
				pos[i] = (int)Math.Round(i / (double)nbins * n);
			}
			double[] result = new double[n];
			lowerQuart = new double[nbins];
			median = new double[nbins];
			upperQuart = new double[nbins];
			binBoundaries = new double[nbins];
			int[] o = ArrayUtil.Order(controlValues);
			int[][] indices = new int[nbins][];
			for (int i = 0; i < nbins; i++) {
				indices[i] = new int[pos[i + 1] - pos[i]];
				Array.Copy(o, pos[i], indices[i], 0, pos[i + 1] - pos[i]);
			}
			for (int i = 0; i < nbins; i++) {
				double[] r = ArrayUtil.SubArray(values, indices[i]);
				double[] c = ArrayUtil.SubArray(controlValues, indices[i]);
				int[] o1 = ArrayUtil.Order(r);
				double rlow = r[o1[(int)Math.Round(0.1587 * (r.Length - 1))]];
				double rmed = r[o1[(int)Math.Round(0.5 * (r.Length - 1))]];
				double rhigh = r[o1[(int)Math.Round(0.8413 * (r.Length - 1))]];
				lowerQuart[i] = rlow;
				median[i] = rmed;
				upperQuart[i] = rhigh;
				binBoundaries[i] = ArrayUtil.Min(c);
				if (indices[i].Length > 2) {
					for (int j = 0; j < indices[i].Length; j++) {
						double ratio = r[j];
						switch (type) {
							case SignificanceType.upper:
								if (rhigh == rmed) {
									result[indices[i][j]] = 1;
								} else {
									double z = (ratio - rmed) / (rhigh - rmed);
									result[indices[i][j]] = Errfunc(z);
								}
								break;
							case SignificanceType.lower:
								if (rlow == rmed) {
									result[indices[i][j]] = 1;
								} else {
									double z = (ratio - rmed) / (rlow - rmed);
									result[indices[i][j]] = Errfunc(z);
								}
								break;
							default:
								if (ratio >= rmed) {
									if (rhigh == rmed) {
										result[indices[i][j]] = 1;
									} else {
										double z = (ratio - rmed) / (rhigh - rmed);
										result[indices[i][j]] = Errfunc(z);
									}
								} else {
									if (rlow == rmed) {
										result[indices[i][j]] = 1;
									} else {
										double z = (ratio - rmed) / (rlow - rmed);
										result[indices[i][j]] = Errfunc(z);
									}
								}
								break;
						}
					}
				} else {
					for (int j = 0; j < indices[i].Length; j++) {
						result[indices[i][j]] = 1;
					}
				}
			}
			return result;
		}

		public static double Errfunc(double z) {
			try {
				return (NumericalRecipes.Erffc(z / Math.Sqrt(2.0))) / 2.0;
			} catch (Exception) {
				return 1;
			}
		}

		public static void CountOutliers(double[] x, out int nOutliers) {
			int n = x.Length;
			if (n < 2) {
				nOutliers = 0;
				return;
			}
			int[] o = ArrayUtil.Order(x);
			int min = (n - 2) / 4;
			double minx;
			if ((n - 2) % 4 == 0) {
				minx = x[o[min]];
			} else {
				minx = 0.5 * (x[o[min]] + x[o[min + 1]]);
			}
			int max = (3 * n - 2) / 4;
			double maxx;
			if ((3 * n - 2) % 4 == 0) {
				maxx = x[o[max]];
			} else {
				maxx = 0.5 * (x[o[max]] + x[o[max + 1]]);
			}
			double upperWhisker = maxx + (maxx - minx) * 1.5;
			double lowerWhisker = minx - (maxx - minx) * 1.5;
			nOutliers = 0;
			foreach (double xx in x) {
				if (xx < lowerWhisker || xx > upperWhisker) {
					nOutliers++;
				}
			}
		}

		#region Nested type: TmpPartition

		private class TmpPartition {
			public List<int> partition;
			public int remainder;

			private TmpPartition() {
			}

			internal TmpPartition(int n) {
				remainder = n;
				partition = new List<int>(n);
			}

			internal TmpPartition add(int a) {
				TmpPartition result = new TmpPartition();
				result.remainder = remainder - a;
				result.partition = new List<int>();
				result.partition.AddRange(partition);
				result.partition.Add(a);
				return result;
			}
		}

		#endregion
	}
}