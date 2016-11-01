/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using MaxQuant.Mol;
using MaxQuant.Num;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Tasks {
	public static class T03SilacAssembly {
		public static Dictionary<int, double>[] Vectorize(double[][] masses, float[][][] intens, int[][][] scanIndices,
		                                                  int charge, double[][][] diffDistr, int[] offsets) {
			int n = diffDistr.Length;
			Dictionary<int, double[][]>[] distr = new Dictionary<int, double[][]>[n];
			for (int i = 0; i < n; i++) {
				if (masses[i] == null) {
					continue;
				}
				Dictionary<int, double[][]> q = GetIsotopePatterns(masses[i], intens[i], scanIndices[i]);
				Dictionary<int, double[][]> r = new Dictionary<int, double[][]>();
				foreach (int scanNum in q.Keys) {
					double[][] d = Molecule.Convolute(q[scanNum], diffDistr[i], 0.2, 1);
					r[scanNum] = d;
				}
				distr[i] = r;
			}
			double minMass = Double.MaxValue;
			for (int i = 0; i < n; i++) {
				if (masses[i] == null) {
					continue;
				}
				foreach (double[][] p in distr[i].Values) {
					if (p[0].Length > 0 && p[0][0] < minMass) {
						minMass = p[0][0];
					}
				}
			}
			Dictionary<int, double>[] w = new Dictionary<int, double>[distr.Length];
			for (int i = 0; i < n; i++) {
				if (distr[i] == null) {
					continue;
				}
				w[i] = new Dictionary<int, double>();
				foreach (int scanIndex in distr[i].Keys) {
					double[][] p = distr[i][scanIndex];
					for (int j = 0; j < p[0].Length; j++) {
						int ind = (int) Math.Round(p[0][j] - minMass);
						w[i][ind + 128 * (scanIndex + offsets[i])] = p[1][j];
					}
				}
			}
			return w;
		}

		public static double FitRatio(Dictionary<int, double> v1, Dictionary<int, double> v2) {
			List<double> vec1 = new List<double>();
			List<double> vec2 = new List<double>();
			foreach (int x in v1.Keys) {
				if (v2.ContainsKey(x)) {
					vec1.Add(v1[x]);
					vec2.Add(v2[x]);
				}
			}
			if (vec1.Count < 3) {
				return double.NaN;
			}
			double[] a1 = vec1.ToArray();
			double[] a2 = vec2.ToArray();
			double dummy1;
			double ratio = NumUtil.MedfitOrigin(a1, a2, out dummy1);
			if (ratio < 0) {
				return double.NaN;
			}
			return ratio;
		}

		public static double Correlate(Dictionary<int, double> v1, IDictionary<int, double> v2) {
			List<double> vec1 = new List<double>();
			List<double> vec2 = new List<double>();
			foreach (int x in v1.Keys) {
				if (v2.ContainsKey(x)) {
					vec1.Add(v1[x]);
					vec2.Add(v2[x]);
				}
			}
			return ArrayUtil.Correlation(vec1.ToArray(), vec2.ToArray());
		}

		public static double CalcRatio(double[] masses1, double[] weights1, double[] masses2, double[] weights2) {
			double minMass = Math.Min(masses1[0], masses2[0]);
			double maxMass = Math.Max(masses1[masses1.Length - 1], masses2[masses2.Length - 1]);
			int len = (int) Math.Round(maxMass - minMass) + 1;
			double[] w1 = new double[len];
			double[] w2 = new double[len];
			for (int j = 0; j < masses1.Length; j++) {
				int ind = (int) Math.Round(masses1[j] - minMass);
				w1[ind] = weights1[j];
			}
			for (int j = 0; j < masses2.Length; j++) {
				int ind = (int) Math.Round(masses2[j] - minMass);
				w2[ind] = weights2[j];
			}
			double dummy1;
			double dummy2;
			double dummy3;
			double ratio;
			NumUtil.FitOrigin(w1, w2, null, out ratio, out dummy1, out dummy2, out dummy3);
			return ratio;
		}

		public static bool IsotopePatternsMatch(double[] masses1, double[] weights1, double[] errors1, double[] masses2,
		                                        double[] weights2, double[] errors2, double sigmas,
		                                        double correlationThreshold, out double isotopeCorrelation) {
			isotopeCorrelation = Double.NaN;
			double minMass1 = masses1[0];
			double minMass2 = masses2[0];
			double maxMass1 = masses1[masses1.Length - 1];
			double maxMass2 = masses2[masses2.Length - 1];
			if (minMass1 > maxMass2 + 0.1) {
				return false;
			}
			if (minMass2 > maxMass1 + 0.1) {
				return false;
			}
			int maxind1 = ArrayUtil.MaxInd(weights1);
			int maxind2 = ArrayUtil.MaxInd(weights2);
			bool m1 = MassMatch(masses1[maxind1], errors1[maxind1], masses2, errors2, sigmas);
			bool m2 = MassMatch(masses2[maxind2], errors2[maxind2], masses1, errors1, sigmas);
			if (!m1 && !m2) {
				return false;
			}
			double minMass = Math.Min(minMass1, minMass2);
			for (int i = 0; i < masses1.Length; i++) {
				masses1[i] -= minMass;
			}
			for (int i = 0; i < masses2.Length; i++) {
				masses2[i] -= minMass;
			}
			int[] intMasses1 = new int[masses1.Length];
			for (int i = 0; i < masses1.Length; i++) {
				intMasses1[i] = (int) Math.Round(masses1[i]);
			}
			int[] intMasses2 = new int[masses2.Length];
			for (int i = 0; i < masses2.Length; i++) {
				intMasses2[i] = (int) Math.Round(masses2[i]);
			}
			int intMassCount = Math.Max(intMasses1[intMasses1.Length - 1],
			                            intMasses2[intMasses2.Length - 1]) + 1;
			double[] profile1 = new double[intMassCount + 2];
			for (int i = 0; i < weights1.Length; i++) {
				profile1[intMasses1[i] + 1] = weights1[i];
			}
			double[] profile2 = new double[intMassCount + 2];
			for (int i = 0; i < weights2.Length; i++) {
				profile2[intMasses2[i] + 1] = weights2[i];
			}
			isotopeCorrelation = ArrayUtil.Cosine(profile1, profile2);
			return isotopeCorrelation >= correlationThreshold;
		}

		public static void Interpolate(ref int[] indices, ref double[] values) {
			int len = indices[indices.Length - 1] - indices[0] + 1;
			int[] newIndices = new int[len];
			double[] newValues = new double[len];
			for (int i = 0; i < len; i++) {
				newIndices[i] = indices[0] + i;
			}
			for (int i = 0; i < indices.Length; i++) {
				int x1 = indices[i] - indices[0];
				newValues[x1] = values[i];
				if (i == indices.Length - 1) {
					continue;
				}
				int x2 = indices[i + 1] - indices[0];
				double v2 = values[i + 1];
				double v1 = values[i];
				double a = (v2 - v1) / (x2 - x1);
				for (int xx = x1 + 1; xx < x2; xx++) {
					newValues[xx] = a * (xx - x1) + v1;
				}
			}
			indices = newIndices;
			values = newValues;
		}

		public static double Correlate(double[] profile1, int start1, double[] profile2, int start2, int tol, out int offset) {
			tol = Math.Abs(tol);
			offset = 0;
			double bestCorr = -double.MaxValue;
			for (int i = -tol; i <= tol; i++) {
				double corr = Correlate(profile1, start1, profile2, start2 + i);
				if (corr > bestCorr) {
					bestCorr = corr;
					offset = i;
				}
			}
			return bestCorr;
		}

		public static SilacPair[,] FindTwoSidedTriangles(ref SilacPair[] pairs01, ref SilacPair[] pairs02,
		                                                 ref SilacPair[] pairs12, HashSet<int> isotopeIndices) {
			List<int[]> candidates = new List<int[]>();
			for (int i = 0; i < pairs01.Length; i++) {
				for (int j = 0; j < pairs02.Length; j++) {
					int ind0 = pairs01[i].Index1;
					if (ind0 == pairs02[j].Index1) {
						int ind1 = pairs01[i].Index2;
						int ind2 = pairs02[j].Index2;
						if (isotopeIndices.Contains(ind0) || isotopeIndices.Contains(ind1) || isotopeIndices.Contains(ind2)) {
							continue;
						}
						if (!pairs01[i].EqualLabeledAaCounts(pairs02[j])) {
							continue;
						}
						candidates.Add(new int[] {i, j, -1});
					}
				}
			}
			for (int i = 0; i < pairs01.Length; i++) {
				for (int j = 0; j < pairs12.Length; j++) {
					int ind1 = pairs01[i].Index2;
					if (ind1 == pairs12[j].Index1) {
						int ind0 = pairs01[i].Index1;
						int ind2 = pairs12[j].Index2;
						if (isotopeIndices.Contains(ind0) || isotopeIndices.Contains(ind1) || isotopeIndices.Contains(ind2)) {
							continue;
						}
						if (!pairs01[i].EqualLabeledAaCounts(pairs12[j])) {
							continue;
						}
						candidates.Add(new int[] {i, -1, j});
					}
				}
			}
			for (int i = 0; i < pairs02.Length; i++) {
				for (int j = 0; j < pairs12.Length; j++) {
					int ind2 = pairs02[i].Index2;
					if (ind2 == pairs12[j].Index2) {
						int ind0 = pairs02[i].Index1;
						int ind1 = pairs12[j].Index1;
						if (isotopeIndices.Contains(ind0) || isotopeIndices.Contains(ind1) || isotopeIndices.Contains(ind2)) {
							continue;
						}
						if (!pairs02[i].EqualLabeledAaCounts(pairs12[j])) {
							continue;
						}
						candidates.Add(new int[] {-1, i, j});
					}
				}
			}
			candidates = MakeConsistent(candidates, pairs01, pairs02, pairs12);
			foreach (int[] c in candidates) {
				int ind0 = -1;
				int ind1 = -1;
				int ind2 = -1;
				if (c[0] >= 0) {
					ind0 = pairs01[c[0]].Index1;
					ind1 = pairs01[c[0]].Index2;
				}
				if (c[1] >= 0) {
					ind0 = pairs02[c[1]].Index1;
					ind2 = pairs02[c[1]].Index2;
				}
				if (c[2] >= 0) {
					ind1 = pairs12[c[2]].Index1;
					ind2 = pairs12[c[2]].Index2;
				}
				isotopeIndices.Add(ind0);
				isotopeIndices.Add(ind1);
				isotopeIndices.Add(ind2);
			}
			SilacPair[,] result = new SilacPair[candidates.Count,3];
			for (int i = 0; i < candidates.Count; i++) {
				int[] c = candidates[i];
				result[i, 0] = (c[0] == -1) ? null : pairs01[c[0]];
				result[i, 1] = (c[1] == -1) ? null : pairs02[c[1]];
				result[i, 2] = (c[2] == -1) ? null : pairs12[c[2]];
			}
			List<SilacPair> newPairs01 = new List<SilacPair>();
			foreach (SilacPair p in pairs01) {
				if (isotopeIndices.Contains(p.Index1) || isotopeIndices.Contains(p.Index2)) {
					continue;
				}
				newPairs01.Add(p);
			}
			pairs01 = newPairs01.ToArray();
			List<SilacPair> newPairs02 = new List<SilacPair>();
			foreach (SilacPair p in pairs02) {
				if (isotopeIndices.Contains(p.Index1) || isotopeIndices.Contains(p.Index2)) {
					continue;
				}
				newPairs02.Add(p);
			}
			pairs02 = newPairs02.ToArray();
			List<SilacPair> newPairs12 = new List<SilacPair>();
			foreach (SilacPair p in pairs12) {
				if (isotopeIndices.Contains(p.Index1) || isotopeIndices.Contains(p.Index2)) {
					continue;
				}
				newPairs12.Add(p);
			}
			pairs12 = newPairs12.ToArray();
			return result;
		}

		public static SilacPair[,] FindFullTriangles(ref SilacPair[] pairs01, ref SilacPair[] pairs02,
		                                             ref SilacPair[] pairs12, HashSet<int> isotopeIndices) {
			List<int[]> firstCandidates = new List<int[]>();
			for (int i = 0; i < pairs01.Length; i++) {
				for (int j = 0; j < pairs02.Length; j++) {
					if (pairs01[i].Index1 == pairs02[j].Index1) {
						if (pairs01[i].EqualLabeledAaCounts(pairs02[j])) {
							firstCandidates.Add(new int[] {i, j});
						}
					}
				}
			}
			List<int[]> candidates = new List<int[]>();
			for (int i = 0; i < firstCandidates.Count; i++) {
				int[] candidate = firstCandidates[i];
				for (int j = 0; j < pairs12.Length; j++) {
					if (pairs01[candidate[0]].Index2 != pairs12[j].Index1) {
						continue;
					}
					if (!pairs01[candidate[0]].EqualLabeledAaCounts(pairs12[j])) {
						continue;
					}
					if (pairs02[candidate[1]].Index2 != pairs12[j].Index2) {
						continue;
					}
					if (!pairs02[candidate[1]].EqualLabeledAaCounts(pairs12[j])) {
						continue;
					}
					int i0 = pairs01[candidate[0]].Index1;
					int i1 = pairs01[candidate[0]].Index2;
					int i2 = pairs02[candidate[1]].Index2;
					if (i0 == i1 || i0 == i2 || i1 == i2) {
						continue;
					}
					candidates.Add(new int[] {candidate[0], candidate[1], j});
				}
			}
			candidates = MakeConsistent(candidates, pairs01, pairs02, pairs12);
			foreach (int[] c in candidates) {
				isotopeIndices.Add(pairs01[c[0]].Index1);
				isotopeIndices.Add(pairs01[c[0]].Index2);
				isotopeIndices.Add(pairs02[c[1]].Index2);
			}
			SilacPair[,] result = new SilacPair[candidates.Count,3];
			for (int i = 0; i < candidates.Count; i++) {
				result[i, 0] = pairs01[candidates[i][0]];
				result[i, 1] = pairs02[candidates[i][1]];
				result[i, 2] = pairs12[candidates[i][2]];
			}
			List<SilacPair> newPairs01 = new List<SilacPair>();
			foreach (SilacPair p in pairs01) {
				if (isotopeIndices.Contains(p.Index1) || isotopeIndices.Contains(p.Index2)) {
					continue;
				}
				newPairs01.Add(p);
			}
			pairs01 = newPairs01.ToArray();
			List<SilacPair> newPairs02 = new List<SilacPair>();
			foreach (SilacPair p in pairs02) {
				if (isotopeIndices.Contains(p.Index1) || isotopeIndices.Contains(p.Index2)) {
					continue;
				}
				newPairs02.Add(p);
			}
			pairs02 = newPairs02.ToArray();
			List<SilacPair> newPairs12 = new List<SilacPair>();
			foreach (SilacPair p in pairs12) {
				if (isotopeIndices.Contains(p.Index1) || isotopeIndices.Contains(p.Index2)) {
					continue;
				}
				newPairs12.Add(p);
			}
			pairs12 = newPairs12.ToArray();
			return result;
		}

		public static double[][] GetIsotopePattern(IsotopeCluster c, out double[] massErrors, double[] centerMz,
		                                           float[] intensities, float[] centerMzErrors) {
			int[] members = c.Members;
			int charge = c.Charge;
			double[][] result = new double[2][];
			result[0] = new double[members.Length];
			result[1] = new double[members.Length];
			massErrors = new double[members.Length];
			for (int i = 0; i < members.Length; i++) {
				result[0][i] = charge * (centerMz[members[i]] - MolUtil.MassProton);
				result[1][i] = intensities[members[i]];
				massErrors[i] = centerMzErrors[members[i]] * charge;
			}
			return result;
		}

		public static double[] SumWeights(double[][][] distr, out int[] offsets) {
			double minMass = Double.MaxValue;
			double maxMass = -Double.MaxValue;
			for (int i = 0; i < distr.Length; i++) {
				if (distr[i] != null && distr[i][0][0] < minMass) {
					minMass = distr[i][0][0];
				}
				if (distr[i] != null && distr[i][0][distr[i][0].Length - 1] > maxMass) {
					maxMass = distr[i][0][distr[i][0].Length - 1];
				}
			}
			int len = (int) Math.Round(maxMass - minMass) + 1;
			double[] weights = new double[len];
			offsets = new int[distr.Length];
			for (int i = 0; i < distr.Length; i++) {
				if (distr[i] == null) {
					continue;
				}
				for (int j = 0; j < distr[i][0].Length; j++) {
					int ind = (int) Math.Round(distr[i][0][j] - minMass);
					if (j == 0) {
						offsets[i] = ind;
					}
					weights[ind] += distr[i][1][j];
				}
			}
			return weights;
		}

		public static SilacPair[] MakeConsistent(SilacPair[] pairs) {
			double[] absLogRatio = new double[pairs.Length];
			for (int i = 0; i < pairs.Length; i++) {
				double r = pairs[i].Ratio;
				absLogRatio[i] = r > 0 ? Math.Abs(Math.Log(r)) : Double.MaxValue;
			}
			int[] o = ArrayUtil.Order(absLogRatio);
			HashSet<int> indices = new HashSet<int>();
			List<SilacPair> result = new List<SilacPair>();
			for (int i = 0; i < pairs.Length; i++) {
				SilacPair p = pairs[o[i]];
				int ind1 = p.Index1;
				int ind2 = p.Index2;
				if ((!indices.Contains(ind1)) && (!indices.Contains(ind2))) {
					result.Add(p);
					indices.Add(ind1);
					indices.Add(ind2);
				}
			}
			return result.ToArray();
		}

		public static int Fit(double[] theorWeights, double[] weights, out double maxCorr) {
			maxCorr = -Double.MaxValue;
			int maxCorrInd = Int32.MinValue;
			for (int i = -theorWeights.Length + 1; i < weights.Length - 1; i++) {
				int start = Math.Min(i, 0);
				int end = Math.Max(i + theorWeights.Length, weights.Length);
				int len = end - start;
				double[] p1 = new double[len];
				double[] p2 = new double[len];
				for (int j = 0; j < theorWeights.Length; j++) {
					p1[i + j - start] = theorWeights[j];
				}
				for (int j = 0; j < weights.Length; j++) {
					p2[j - start] = weights[j];
				}
				double corr = ArrayUtil.Cosine(p1, p2);
				if (corr > maxCorr) {
					maxCorr = corr;
					maxCorrInd = i;
				}
			}
			return -maxCorrInd;
		}

		public static double GetMassErrorEstimate(int[] isoClusterIndices, int[] isotopePatternStart,
		                                          double[] massShifts, IsotopeCluster[] isotopeClusters, IPeakList peakList) {
			double[] m = new double[ArrayUtil.nBoots];
			for (int i = 0; i < ArrayUtil.nBoots; i++) {
				m[i] = GetFullIsotopePatternMassEstimateBootstrap(isoClusterIndices, isotopePatternStart, massShifts,
				                                                  isotopeClusters, peakList);
			}
			double mean = 0;
			for (int i = 0; i < ArrayUtil.nBoots; i++) {
				mean += m[i];
			}
			mean /= ArrayUtil.nBoots;
			double err = 0;
			for (int i = 0; i < ArrayUtil.nBoots; i++) {
				err += (mean - m[i]) * (mean - m[i]);
			}
			err /= ArrayUtil.nBoots;
			err = Math.Sqrt(err);
			return err;
		}

		private static Dictionary<int, double>[] Vectorize(IsotopeCluster[] ics, double[][][] diffDistr, int[] offsets,
		                                                   double[] centerMz, IPeakList peakList) {
			int n = diffDistr.Length;
			double[][] masses = new double[n][];
			float[][][] intens = new float[n][][];
			int[][][] scanIndices = new int[n][][];
			int charge = -1;
			for (int i = 0; i < n; i++) {
				if (ics[i] == null) {
					continue;
				}
				int[] members = ics[i].Members;
				charge = ics[i].Charge;
				masses[i] = new double[members.Length];
				intens[i] = new float[members.Length][];
				scanIndices[i] = new int[members.Length][];
				for (int ij = 0; ij < members.Length; ij++) {
					masses[i][ij] = charge * (centerMz[members[ij]] - MolUtil.MassProton);
					Peak p = peakList.GetPeakKeep(members[ij]);
					intens[i][ij] = p.GetOriginalIntensities();
					scanIndices[i][ij] = p.GetScanIndices();
				}
			}
			return Vectorize(masses, intens, scanIndices, charge, diffDistr, offsets);
		}

		private static Dictionary<int, double[][]> GetIsotopePatterns(double[] masses, float[][] intens, int[][] scanIndices) {
			HashSet<int> nums = new HashSet<int>();
			for (int i = 0; i < scanIndices.Length; i++) {
				foreach (int x in scanIndices[i]) {
					nums.Add(x);
				}
			}
			int[] allScanNums = nums.ToArray();
			Array.Sort(allScanNums);
			Dictionary<int, double[][]> profiles = new Dictionary<int, double[][]>();
			for (int i = 0; i < allScanNums.Length; i++) {
				int scanNum = allScanNums[i];
				List<double> xx = new List<double>();
				List<double> yy = new List<double>();
				for (int j = 0; j < scanIndices.Length; j++) {
					int ind = Array.BinarySearch(scanIndices[j], scanNum);
					if (ind >= 0) {
						xx.Add(masses[j]);
						yy.Add(intens[j][ind]);
					}
				}
				double[][] profile = new double[2][];
				profile[0] = xx.ToArray();
				profile[1] = yy.ToArray();
				profiles[scanNum] = profile;
			}
			return profiles;
		}

		private static Dictionary<int, double>[] Vectorize(int[] isoInd, double[][][] diffDistr, int[] offsets,
		                                                   IPeakList peakList, IsotopeCluster[] isotopeClusters,
		                                                   double[] centerMz) {
			IsotopeCluster[] ics = new IsotopeCluster[isoInd.Length];
			for (int i = 0; i < ics.Length; i++) {
				if (isoInd[i] != -1) {
					ics[i] = isotopeClusters[isoInd[i]];
				}
			}
			return Vectorize(ics, diffDistr, offsets, centerMz, peakList);
		}

		private static double Correlate(int[] isoInd, double[][][] diffDistr, int[] offsets, IPeakList peakList,
		                                IsotopeCluster[] isotopeClusters, double[] centerMz) {
			Dictionary<int, double>[] w = Vectorize(isoInd, diffDistr, offsets, peakList, isotopeClusters, centerMz);
			return Correlate(w[0], w[1]);
		}

		private static bool MassMatch(double m1, double e1, double m2, double e2, double sigmas) {
			double diff = m1 - m2;
			return diff * diff <= Math.Max(sigmas * sigmas * (e1 * e1 + e2 * e2), 0.002 * 0.002);
		}

		private static bool MassMatch(double mass, double error, double[] masses, double[] errors, double sigmas) {
			for (int i = 0; i < masses.Length; i++) {
				if (MassMatch(mass, error, masses[i], errors[i], sigmas)) {
					return true;
				}
			}
			return false;
		}

		private static bool IsSilacPair(int isoId1, int isoId2, double[][] pattern1, double[] errors1,
		                                double[][] pattern2, double[] errors2, LabelCombinations labelCombinations,
		                                double isotopeCorrelationThreshold, double sigmas, int offset,
		                                out int[] counts, out double isotopeCorrelation,
		                                out double ratio, IPeakList peakList, IsotopeCluster[] isotopeClusters,
		                                double[] centerMz) {
			double bestCorr = -double.MaxValue;
			int bestCorrInd = -1;
			double[] bestMassesX1 = null;
			double[] bestWeightsX1 = null;
			double[] bestMassesX2 = null;
			double[] bestWeightsX2 = null;
			for (int i = 0; i < labelCombinations.PartitionCount; i++) {
				double[][] d1 = labelCombinations.GetIsotopeDistribution1(i);
				double[][] d2 = labelCombinations.GetIsotopeDistribution2(i);
				double[] massesX1;
				double[] weightsX1;
				double[] errorsX1;
				double[] massesX2;
				double[] weightsX2;
				double[] errorsX2;
				Molecule.ConvoluteWithErrors(pattern1[0], pattern1[1], errors1, d2[0], d2[1], 1.0, 1.0, out massesX1, out weightsX1,
				                             out errorsX1);
				Molecule.ConvoluteWithErrors(pattern2[0], pattern2[1], errors2, d1[0], d1[1], 1.0, 1.0, out massesX2, out weightsX2,
				                             out errorsX2);
				if (
					IsotopePatternsMatch(massesX1, weightsX1, errorsX1, massesX2, weightsX2, errorsX2, sigmas,
					                     isotopeCorrelationThreshold, out isotopeCorrelation)) {
					double corr = Correlate(new int[] {isoId1, isoId2}, new double[][][] {d2, d1}, new int[] {0, offset}, peakList,
					                        isotopeClusters, centerMz);
					if (corr >= isotopeCorrelationThreshold && corr > bestCorr) {
						bestCorr = corr;
						bestCorrInd = i;
						bestMassesX1 = massesX1;
						bestWeightsX1 = weightsX1;
						bestMassesX2 = massesX2;
						bestWeightsX2 = weightsX2;
					}
				}
			}
			if (bestCorrInd != -1) {
				double[][] d1 = labelCombinations.GetIsotopeDistribution1(bestCorrInd);
				double[][] d2 = labelCombinations.GetIsotopeDistribution2(bestCorrInd);
				counts = labelCombinations.GetPartition(bestCorrInd);
				ratio = CalcRatio(bestMassesX1, bestWeightsX1, bestMassesX2, bestWeightsX2);
				isotopeCorrelation = Correlate(new int[] {isoId1, isoId2}, new double[][][] {d2, d1}, new int[] {0, offset},
				                               peakList, isotopeClusters, centerMz);
				return true;
			}
			counts = null;
			isotopeCorrelation = Double.NaN;
			ratio = Double.NaN;
			return false;
		}

		private static double[] GetBoundingBox(IsotopeCluster c, double[] centerMz, float[] minTimes, float[] maxTimes) {
			int[] m = c.Members;
			double minMass = centerMz[m[0]];
			double maxMass = centerMz[m[m.Length - 1]];
			double minRT = Double.MaxValue;
			double maxRT = -Double.MaxValue;
			foreach (int ind in m) {
				if (minTimes[ind] < minRT) {
					minRT = minTimes[ind];
				}
				if (maxTimes[ind] > maxRT) {
					maxRT = maxTimes[ind];
				}
			}
			return new double[] {minMass, maxMass, minRT, maxRT};
		}

		public static double[] CalcIsotopeClusterProfile(IsotopeCluster cl, out int minIndex, IPeakList peakList) {
			int[] m = cl.Members;
			int n = m.Length;
			Peak[] peak = new Peak[n];
			for (int i = 0; i < n; i++) {
				peak[i] = peakList.GetPeakKeep(m[i]);
			}
			int[][] scanIndices = new int[n][];
			int[] firstIndices = new int[n];
			int[] lastIndices = new int[n];
			for (int i = 0; i < n; i++) {
				scanIndices[i] = peak[i].GetScanIndices();
				firstIndices[i] = scanIndices[i][0];
				lastIndices[i] = scanIndices[i][scanIndices[i].Length - 1];
			}
			minIndex = ArrayUtil.Min(firstIndices);
			int maxIndex = ArrayUtil.Max(lastIndices);
			int len = maxIndex - minIndex + 1;
			double[] profile = new double[len + 2];
			for (int i = 0; i < n; i++) {
				int[] indices = scanIndices[i];
				double[] values = new double[indices.Length];
				for (int j = 0; j < values.Length; j++) {
					values[j] = peak[i].GetSmoothIntensity(j);
				}
				Interpolate(ref indices, ref values);
				for (int a = 0; a < scanIndices[i].Length; a++) {
					profile[indices[a] - minIndex + 1] += values[a];
				}
			}
			minIndex -= 1;
			return profile;
		}

		public static double Correlate(double[] profile1, int start1, double[] profile2, int start2) {
			int start = Math.Min(start1, start2);
			int end = Math.Max(start1 + profile1.Length, start2 + profile2.Length);
			int len = end - start;
			double[] p1 = new double[len];
			double[] p2 = new double[len];
			for (int i = 0; i < profile1.Length; i++) {
				p1[start1 - start + i] = profile1[i];
			}
			for (int i = 0; i < profile2.Length; i++) {
				p2[start2 - start + i] = profile2[i];
			}
			return ArrayUtil.Cosine(p1, p2);
		}

		private static SilacPair[] FindPairs(int charge, int minCharge, LabelCombinations labelCombinations, double sigmas,
		                                     double timeCorrelationThreshold, double isotopeCorrelationThreshold,
		                                     int maxTolerance,
		                                     IPeakList peakList, IsotopeCluster[] isotopeClusters, double[] centerMz,
		                                     int[][] isotopeClusterIdsByCharge, float[] minTimes, float[] maxTimes,
		                                     float[] intensities, float[] centerMzErrors) {
			int[] isotopeClustersIds = isotopeClusterIdsByCharge[charge - minCharge];
			double[][][] isotopePatterns = new double[isotopeClustersIds.Length][][];
			double[][] massErrors = new double[isotopeClustersIds.Length][];
			for (int i = 0; i < isotopeClustersIds.Length; i++) {
				isotopePatterns[i] = GetIsotopePattern(isotopeClusters[isotopeClustersIds[i]], out massErrors[i],
				                                       centerMz, intensities, centerMzErrors);
			}
			double[][] boundingBoxes = new double[isotopeClustersIds.Length][];
			for (int i = 0; i < isotopeClustersIds.Length; i++) {
				boundingBoxes[i] = GetBoundingBox(isotopeClusters[isotopeClustersIds[i]], centerMz, minTimes, maxTimes);
			}
			double[][] timeProfiles = new double[isotopeClustersIds.Length][];
			int[] startIndices = new int[isotopeClustersIds.Length];
			bool[] invalid = new bool[isotopeClustersIds.Length];
			for (int i = 0; i < isotopeClustersIds.Length; i++) {
				IsotopeCluster ic = isotopeClusters[isotopeClustersIds[i]];
				timeProfiles[i] = CalcIsotopeClusterProfile(ic, out startIndices[i], peakList);
				invalid[i] = ic.PolymerIndex != -1 || ic.ContaminantIndex != -1;
			}
			double massDiff = (labelCombinations.MaxMassDiff + 1) / charge;
			List<SilacPair> result = new List<SilacPair>();
			for (int i = 0; i < isotopeClustersIds.Length; i++) {
				double minMassI = boundingBoxes[i][0];
				double maxMassI = boundingBoxes[i][1];
				double minRtI = boundingBoxes[i][2];
				double maxRtI = boundingBoxes[i][3];
				for (int j = 0; j < isotopeClustersIds.Length; j++) {
					if (i == j) {
						continue;
					}
					if (invalid[i] || invalid[j]) {
						continue;
					}
					double minMassJ = boundingBoxes[j][0];
					double maxMassJ = boundingBoxes[j][1];
					double minRtJ = boundingBoxes[j][2];
					double maxRtJ = boundingBoxes[j][3];
					if (minRtI > maxRtJ) {
						continue;
					}
					if (minRtJ > maxRtI) {
						continue;
					}
					if (minMassI > maxMassJ + massDiff) {
						continue;
					}
					if (minMassJ > maxMassI + massDiff) {
						continue;
					}
					int offset;
					double timeCorrelation = Correlate(timeProfiles[i], startIndices[i], timeProfiles[j],
					                                   startIndices[j], maxTolerance, out offset);
					if (timeCorrelation < timeCorrelationThreshold) {
						continue;
					}
					int[] counts;
					double isotopeCorrelation;
					double ratio;
					if (IsSilacPair(isotopeClustersIds[i], isotopeClustersIds[j], isotopePatterns[i],
					                massErrors[i], isotopePatterns[j], massErrors[j], labelCombinations,
					                isotopeCorrelationThreshold, sigmas, offset, out counts, out isotopeCorrelation,
					                out ratio, peakList, isotopeClusters, centerMz)) {
						result.Add(new SilacPair(isotopeClustersIds[i], isotopeClustersIds[j], counts,
						                         isotopeCorrelation, timeCorrelation, ratio, offset));
					}
				}
			}
			return result.ToArray();
		}

		private static List<int[]> MakeConsistent(IList<int[]> candidates, SilacPair[] pairs01,
		                                          SilacPair[] pairs02, SilacPair[] pairs12) {
			double[] r = new double[candidates.Count];
			for (int i = 0; i < r.Length; i++) {
				int[] c = candidates[i];
				r[i] = (c[0] == -1 ? 0 : AbsLogRatio(pairs01[c[0]])) +
				       (c[1] == -1 ? 0 : AbsLogRatio(pairs02[c[1]])) +
				       (c[2] == -1 ? 0 : AbsLogRatio(pairs12[c[2]]));
			}
			int[] order = ArrayUtil.Order(r);
			HashSet<int> allInds = new HashSet<int>();
			List<int[]> result = new List<int[]>();
			for (int i = 0; i < order.Length; i++) {
				int[] c = candidates[order[i]];
				int ind0 = -1;
				int ind1 = -1;
				int ind2 = -1;
				if (c[0] >= 0) {
					ind0 = pairs01[c[0]].Index1;
					ind1 = pairs01[c[0]].Index2;
				}
				if (c[1] >= 0) {
					ind0 = pairs02[c[1]].Index1;
					ind2 = pairs02[c[1]].Index2;
				}
				if (c[2] >= 0) {
					ind1 = pairs12[c[2]].Index1;
					ind2 = pairs12[c[2]].Index2;
				}
				if (allInds.Contains(ind0) || allInds.Contains(ind1) || allInds.Contains(ind2)) {
					continue;
				}
				result.Add(c);
				if (ind0 >= 0) {
					allInds.Add(ind0);
				}
				if (ind1 >= 0) {
					allInds.Add(ind1);
				}
				if (ind2 >= 0) {
					allInds.Add(ind2);
				}
			}
			return result;
		}

		private static double AbsLogRatio(SilacPair pair) {
			if (pair.Ratio > 0) {
				return Math.Abs(Math.Log(pair.Ratio));
			}
			return double.MaxValue;
		}

		public static double GetIsotopeClusterIntensity(IsotopeCluster ic, float[] intensities) {
			double result = 0;
			int[] m1 = ic.Members;
			for (int j = 0; j < m1.Length; j++) {
				double x = intensities[m1[j]];
				//if (x > result) {
				//    result = x;
				//}
				result += x;
			}
			return result;
		}

		private static void CalcRatiosDoublets(int[] isoInd, double[][][] diffDistr, int offset, out double ratio10,
		                                       IPeakList peakList, IsotopeCluster[] isotopeClusters,
		                                       double[] centerMz) {
			Dictionary<int, double>[] w = Vectorize(isoInd, diffDistr, new int[] {0, offset}, peakList, isotopeClusters,
			                                        centerMz);
			if (w[0] != null && w[1] != null) {
				ratio10 = FitRatio(w[0], w[1]);
			} else {
				ratio10 = w[0] == null ? double.PositiveInfinity : 0;
			}
		}

		private static double BootstrapWeightedMean(double[] x, double[] weights) {
			int[] ind = ArrayUtil.GetBootstrapIndices(x.Length);
			return WeightedMean(ArrayUtil.SubArray(x, ind), ArrayUtil.SubArray(weights, ind));
		}

		private static double WeightedMean(double[] x, double[] weights) {
			double result = 0;
			double norm = 0;
			for (int i = 0; i < x.Length; i++) {
				result += x[i] * weights[i];
				norm += weights[i];
			}
			return result / norm;
		}

		private static double GetFullIsotopePatternMassEstimateBootstrap(int[] isoClusterIndices, int[] isotopePatternStart,
		                                                                 double[] massShifts, IsotopeCluster[] isotopeClusters,
		                                                                 IPeakList peakList) {
			List<double> masses = new List<double>();
			List<double> weights = new List<double>();
			for (int k = 0; k < isoClusterIndices.Length; k++) {
				if (isoClusterIndices[k] < 0) {
					continue;
				}
				IsotopeCluster ic = isotopeClusters[isoClusterIndices[k]];
				int charge = ic.Charge;
				for (int i = 0; i < ic.Count; i++) {
					if (i + isotopePatternStart[k] >= 0) {
						int ind = ic.Members[i];
						Peak peak = peakList.GetPeakKeep(ind);
						for (int index0 = 0; index0 < peak.Count; index0++) {
							weights.Add(peak.OrigIntensity[index0]);
							double mz = peak.CenterMz[index0];
							double m = (mz - MolUtil.MassProton) * charge;
							masses.Add((m - MolUtil.GetAverageDifferenceToMonoisotope(m, i + isotopePatternStart[k]) - massShifts[k]));
						}
					}
				}
			}
			if (masses.Count > 0) {
				return BootstrapWeightedMean(masses.ToArray(), weights.ToArray());
			}
			for (int k = 0; k < isoClusterIndices.Length; k++) {
				if (isoClusterIndices[k] < 0) {
					continue;
				}
				IsotopeCluster ic = isotopeClusters[isoClusterIndices[k]];
				int charge = ic.Charge;
				for (int i = 0; i < ic.Count; i++) {
					int ind = ic.Members[i];
					Peak peak = peakList.GetPeakKeep(ind);
					for (int index0 = 0; index0 < peak.Count; index0++) {
						weights.Add(peak.OrigIntensity[index0]);
						double mz = peak.CenterMz[index0];
						double m = (mz - MolUtil.MassProton) * charge;
						masses.Add(
							(m - (i + isotopePatternStart[k]) * MolUtil.GetAverageDifferenceToMonoisotope(m, 1) -
							 massShifts[k]));
					}
				}
			}
			return BootstrapWeightedMean(masses.ToArray(), weights.ToArray());
		}

		public static SilacCluster[] FindSilacDoublets(int minCharge, int maxCharge, SilacLabel[] labels1, int maxSilacAa,
		                                               double sigmas, double silacTimeCorrelationThreshold,
		                                               double silacIsotopeCorrelationThreshold, int maxTolerance,
		                                               IPeakList peakList, SilacType silacType,
		                                               IsotopeCluster[] isotopeClusters, double[] centerMz,
		                                               int[][] isotopeClusterIdsByCharge, float[] minTimes, float[] maxTimes,
		                                               float[] intensities, float[] centerMzErrors, AminoAcid[] allAas) {
			List<SilacCluster> silacCl = new List<SilacCluster>();
			int nCharge = maxCharge - minCharge + 1;
			SilacPair[][] pairs1 = new SilacPair[nCharge][];
			char[] allAaLetts = AminoAcid.GetSingleLetters(allAas);
			LabelCombinations labelCombinations1 = new LabelCombinations(labels1, maxSilacAa, allAaLetts);
			for (int i = 0; i < nCharge; i++) {
				pairs1[i] = FindPairs(minCharge + i, minCharge, labelCombinations1, sigmas,
				                      silacTimeCorrelationThreshold, silacIsotopeCorrelationThreshold,
				                      maxTolerance, peakList, isotopeClusters, centerMz, isotopeClusterIdsByCharge,
				                      minTimes, maxTimes, intensities, centerMzErrors);
				pairs1[i] = MakeConsistent(pairs1[i]);
				for (int j = 0; j < pairs1[i].Length; j++) {
					SilacPair pair = pairs1[i][j];
					SilacCluster c = new SilacCluster(pair);
					double[] massErrors;
					int[] offsets;
					double maxCorr;
					double ratio10;
					Molecule diff2 = labelCombinations1.CalcDiff1(c.AaCounts);
					Molecule diff1 = labelCombinations1.CalcDiff2(c.AaCounts);
					double massDiff = diff1.GetMostLikelyMass(0.2) - diff2.GetMostLikelyMass(0.2);
					double[][] p1 = GetIsotopePattern(isotopeClusters[c.IsotopeClusterindex1], out massErrors, centerMz,
					                                  intensities, centerMzErrors);
					double[][] p2 = GetIsotopePattern(isotopeClusters[c.IsotopeClusterindex2], out massErrors, centerMz,
					                                  intensities, centerMzErrors);
					double[][] d1 = diff1.GetIsotopeDistribution(0.2);
					double[][] d2 = diff2.GetIsotopeDistribution(0.2);
					double approxMass = p1[0][0];
					double[][] isotopePattern1 = Molecule.Convolute(p1, d1, 0.2, 1);
					double[][] isotopePattern2 = Molecule.Convolute(p2, d2, 0.2, 1);
					double[] weights = SumWeights(new double[][][] {isotopePattern1, isotopePattern2}, out offsets);
					CalcRatiosDoublets(c.GetIsotopeClusterIndices(silacType), new double[][][] {d1, d2},
					                   pair.Offset, out ratio10, peakList, isotopeClusters, centerMz);
					double[][] averagePattern = MolUtil.GetAverageIsotopePattern(approxMass, true);
					averagePattern = Molecule.Convolute(averagePattern, d1, 0.2, 1e-6);
					int maxCorrInd = Fit(averagePattern[1], weights, out maxCorr);
					c.TheorIsotopeCorr = maxCorr;
					c.IsotopePatternStart1 = maxCorrInd + offsets[0];
					c.IsotopePatternStart2 = maxCorrInd + offsets[1];
					c.Ratio10 = ratio10;
					c.Intensity1 = GetIsotopeClusterIntensity(isotopeClusters[c.IsotopeClusterindex1], intensities);
					c.Intensity2 = GetIsotopeClusterIntensity(isotopeClusters[c.IsotopeClusterindex2], intensities);
					c.MassDiff1 = massDiff;
					c.SetMassError(GetMassErrorEstimate(c.GetIsotopeClusterIndices(silacType),
					                                    c.GetIsotopePatternStarts(silacType),
					                                    c.GetMassDiffs(silacType), isotopeClusters, peakList));
					silacCl.Add(c);
				}
			}
			return silacCl.ToArray();
		}

		public static SilacCluster[] FindSilacTriplets(int minCharge, int maxCharge, SilacLabel[] labels1,
		                                               SilacLabel[] labels2, int maxSilacAa, double sigmas,
		                                               double silacTimeCorrelationThreshold,
		                                               double silacIsotopeCorrelationThreshold, int maxTolerance,
		                                               IPeakList peakList, SilacType silacType,
		                                               IsotopeCluster[] isotopeClusters, double[] centerMz,
		                                               int[][] isotopeClusterIdsByCharge, float[] minTimes, float[] maxTimes,
		                                               float[] intensities, float[] centerMzErrors, AminoAcid[] allAas) {
			int nCharge = maxCharge - minCharge + 1;
			char[] allAaLetts = AminoAcid.GetSingleLetters(allAas);
			LabelCombinations labelCombinations01 = new LabelCombinations(labels1, maxSilacAa, allAaLetts);
			LabelCombinations labelCombinations02 = new LabelCombinations(labels2, maxSilacAa, allAaLetts);
			LabelCombinations labelCombinations12 =
				new LabelCombinations(labels1, labels2, maxSilacAa, allAaLetts);
			List<SilacCluster> allClusters = new List<SilacCluster>();
			for (int i = 0; i < nCharge; i++) {
				HashSet<int> isotopeIndices = new HashSet<int>();
				SilacPair[] pairs01 = FindPairs(minCharge + i, minCharge, labelCombinations01,
				                                sigmas, silacTimeCorrelationThreshold,
				                                silacIsotopeCorrelationThreshold,
				                                maxTolerance, peakList, isotopeClusters, centerMz,
				                                isotopeClusterIdsByCharge, minTimes,
				                                maxTimes, intensities, centerMzErrors);
				SilacPair[] pairs02 = FindPairs(minCharge + i, minCharge, labelCombinations02,
				                                sigmas, silacTimeCorrelationThreshold,
				                                silacIsotopeCorrelationThreshold,
				                                maxTolerance, peakList, isotopeClusters, centerMz,
				                                isotopeClusterIdsByCharge, minTimes,
				                                maxTimes, intensities, centerMzErrors);
				SilacPair[] pairs12 = FindPairs(minCharge + i, minCharge, labelCombinations12,
				                                sigmas, silacTimeCorrelationThreshold,
				                                silacIsotopeCorrelationThreshold,
				                                maxTolerance, peakList, isotopeClusters, centerMz,
				                                isotopeClusterIdsByCharge, minTimes,
				                                maxTimes, intensities, centerMzErrors);
				SilacPair[,] fullTriangles =
					FindFullTriangles(ref pairs01, ref pairs02, ref pairs12, isotopeIndices);
				for (int j = 0; j < fullTriangles.GetLength(0); j++) {
					allClusters.Add(new SilacCluster(fullTriangles[j, 0], fullTriangles[j, 1], fullTriangles[j, 2]));
				}
				SilacPair[,] twoSidedTriangles =
					FindTwoSidedTriangles(ref pairs01, ref pairs02, ref pairs12, isotopeIndices);
				for (int j = 0; j < twoSidedTriangles.GetLength(0); j++) {
					allClusters.Add(new SilacCluster(twoSidedTriangles[j, 0], twoSidedTriangles[j, 1], twoSidedTriangles[j, 2]));
				}
			}
			foreach (SilacCluster c in allClusters) {
				CalcSilacTripletProperties(c, labelCombinations01, labelCombinations02, isotopeClusters, centerMz, centerMzErrors,
				                           intensities, silacType, peakList);
			}
			return allClusters.ToArray();
		}

		private static void CalcSilacTripletProperties(SilacCluster c, LabelCombinations labelCombinations01,
		                                               LabelCombinations labelCombinations02, IsotopeCluster[] isotopeClusters,
		                                               double[] centerMz, float[] centerMzErrors, float[] intensities,
		                                               SilacType silacType, IPeakList peakList) {
			double[] massErrors;
			int[] offsets;
			double maxCorr;
			double ratio10;
			double ratio20;
			double ratio21;
			Molecule x1 = labelCombinations01.CalcDiff1(c.AaCounts);
			Molecule x2 = labelCombinations01.CalcDiff2(c.AaCounts);
			double massDiff1 = x2.GetMostLikelyMass(0.2) - x1.GetMostLikelyMass(0.2);
			Molecule y1 = labelCombinations02.CalcDiff1(c.AaCounts);
			Molecule y2 = labelCombinations02.CalcDiff2(c.AaCounts);
			double massDiff2 = y2.GetMostLikelyMass(0.2) - y1.GetMostLikelyMass(0.2);
			Molecule max1 = Molecule.Max(x1, y1);
			Molecule x11 = Molecule.GetDifferences(max1, x1)[0];
			Molecule y11 = Molecule.GetDifferences(max1, y1)[0];
			Molecule diff1 = max1;
			Molecule diff2 = Molecule.Sum(x2, x11);
			Molecule diff3 = Molecule.Sum(y2, y11);
			Molecule max = Molecule.Max(Molecule.Max(diff1, diff2), diff3);
			diff1 = Molecule.GetDifferences(max, diff1)[0];
			diff2 = Molecule.GetDifferences(max, diff2)[0];
			diff3 = Molecule.GetDifferences(max, diff3)[0];
			double[][] p1 = null;
			double[][] p2 = null;
			double[][] p3 = null;
			if (c.IsotopeClusterindex1 >= 0) {
				p1 = GetIsotopePattern(isotopeClusters[c.IsotopeClusterindex1], out massErrors, centerMz, intensities,
				                       centerMzErrors);
			}
			if (c.IsotopeClusterindex2 >= 0) {
				p2 = GetIsotopePattern(isotopeClusters[c.IsotopeClusterindex2], out massErrors, centerMz, intensities,
				                       centerMzErrors);
			}
			if (c.IsotopeClusterindex3 >= 0) {
				p3 = GetIsotopePattern(isotopeClusters[c.IsotopeClusterindex3], out massErrors, centerMz, intensities,
				                       centerMzErrors);
			}
			double[][] d1 = diff1.GetIsotopeDistribution(0.2);
			double[][] d2 = diff2.GetIsotopeDistribution(0.2);
			double[][] d3 = diff3.GetIsotopeDistribution(0.2);
			double[][] isotopePattern1 = null;
			double[][] isotopePattern2 = null;
			double[][] isotopePattern3 = null;
			if (p1 != null) {
				isotopePattern1 = Molecule.Convolute(p1, d1, 0.2, 1);
			}
			if (p2 != null) {
				isotopePattern2 = Molecule.Convolute(p2, d2, 0.2, 1);
			}
			if (p3 != null) {
				isotopePattern3 = Molecule.Convolute(p3, d3, 0.2, 1);
			}
			double[] weights =
				SumWeights(new double[][][] {isotopePattern1, isotopePattern2, isotopePattern3}, out offsets);
			CalcRatiosTriplets(c.GetIsotopeClusterIndices(silacType), new double[][][] {d1, d2, d3},
			                   c.GetOffsets(silacType), out ratio10, out ratio20, out ratio21, peakList, isotopeClusters,
			                   centerMz);
			double approxMass = double.NaN;
			if (p1 != null) {
				approxMass = p1[0][0];
			} else if (p2 != null) {
				approxMass = p2[0][0];
			} else if (p3 != null) {
				approxMass = p3[0][0];
			}
			double[][] averagePattern = MolUtil.GetAverageIsotopePattern(approxMass, true);
			averagePattern = Molecule.Convolute(averagePattern, d1, 0.2, 1e-6);
			int maxCorrInd = Fit(averagePattern[1], weights, out maxCorr);
			c.TheorIsotopeCorr = maxCorr;
			if (c.IsotopeClusterindex1 >= 0) {
				c.IsotopePatternStart1 = maxCorrInd + offsets[0];
			} else {
				c.IsotopePatternStart1 = int.MaxValue;
			}
			if (c.IsotopeClusterindex2 >= 0) {
				c.IsotopePatternStart2 = maxCorrInd + offsets[1];
			} else {
				c.IsotopePatternStart2 = int.MaxValue;
			}
			if (c.IsotopeClusterindex3 >= 0) {
				c.IsotopePatternStart3 = maxCorrInd + offsets[2];
			} else {
				c.IsotopePatternStart3 = int.MaxValue;
			}
			c.Ratio10 = ratio10;
			c.Ratio20 = ratio20;
			c.Ratio21 = ratio21;
			c.Intensity1 = GetIsotopeClusterIntensity(isotopeClusters[c.IsotopeClusterindex1], intensities);
			c.Intensity2 = GetIsotopeClusterIntensity(isotopeClusters[c.IsotopeClusterindex2], intensities);
			c.Intensity3 = GetIsotopeClusterIntensity(isotopeClusters[c.IsotopeClusterindex3], intensities);
			c.MassDiff1 = massDiff1;
			c.MassDiff2 = massDiff2;
			c.SetMassError(GetMassErrorEstimate(c.GetIsotopeClusterIndices(silacType),
			                                    c.GetIsotopePatternStarts(silacType),
			                                    c.GetMassDiffs(silacType), isotopeClusters, peakList));
		}

		private static void CalcRatiosTriplets(int[] isoInd, double[][][] diffDistr, int[] offsets, out double ratio10,
		                                       out double ratio20, out double ratio21, IPeakList peakList,
		                                       IsotopeCluster[] isotopeClusters, double[] centerMz) {
			if (isoInd[0] != -1 && isoInd[1] != -1) {
				Dictionary<int, double>[] w =
					Vectorize(new int[] {isoInd[0], isoInd[1]}, new double[][][] {diffDistr[0], diffDistr[1]},
					          new int[] {0, offsets[0]}, peakList, isotopeClusters, centerMz);
				ratio10 = FitRatio(w[0], w[1]);
			} else {
				ratio10 = isoInd[0] == -1 ? double.PositiveInfinity : 0;
			}
			if (isoInd[0] != -1 && isoInd[2] != -1) {
				Dictionary<int, double>[] w =
					Vectorize(new int[] {isoInd[0], isoInd[2]}, new double[][][] {diffDistr[0], diffDistr[2]},
					          new int[] {0, offsets[1]}, peakList, isotopeClusters, centerMz);
				ratio20 = FitRatio(w[0], w[1]);
			} else {
				ratio20 = isoInd[0] == -1 ? double.PositiveInfinity : 0;
			}
			if (isoInd[1] != -1 && isoInd[2] != -1) {
				Dictionary<int, double>[] w =
					Vectorize(new int[] {isoInd[1], isoInd[2]}, new double[][][] {diffDistr[1], diffDistr[2]},
					          new int[] {0, offsets[2]}, peakList, isotopeClusters, centerMz);
				ratio21 = FitRatio(w[0], w[1]);
			} else {
				ratio21 = isoInd[1] == -1 ? double.PositiveInfinity : 0;
			}
		}
	}
}