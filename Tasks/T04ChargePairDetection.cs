/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Collections.Generic;
using MaxQuant.Mol;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Tasks {
	public static class T04ChargePairDetection {
		public static double[] CalcSilacClusterProfile(SilacCluster cl, out int minIndex, SilacType type,
		                                               IsotopeCluster[] isotopeClusters, IPeakList peakList,
		                                               float[] intensities) {
			int[] m = cl.GetIsotopeClusterIndices(type);
			int n = m.Length;
			if (n == 1) {
				return T03SilacAssembly.CalcIsotopeClusterProfile(isotopeClusters[m[0]], out minIndex, peakList);
			}
			double[][] profiles = new double[n][];
			int[] mins = new int[n];
			int[] ends = new int[n];
			double[] intens = new double[n];
			for (int i = 0; i < n; i++) {
				if (m[i] == -1) {
					continue;
				}
				IsotopeCluster ic = isotopeClusters[m[i]];
				profiles[i] = T03SilacAssembly.CalcIsotopeClusterProfile(ic, out mins[i], peakList);
				intens[i] = T03SilacAssembly.GetIsotopeClusterIntensity(ic, intensities);
				ends[i] = mins[i] + profiles[i].Length;
			}
			minIndex = ArrayUtil.Min(mins);
			int maxEnd = ArrayUtil.Max(ends);
			int len = maxEnd - minIndex;
			double[] result = new double[len];
			double norm = 0;
			for (int i = 0; i < n; i++) {
				if (m[i] == -1) {
					continue;
				}
				norm += intens[i];
				for (int j = 0; j < profiles[i].Length; j++) {
					result[j + mins[i] - minIndex] += profiles[i][j] * intens[i];
				}
			}
			for (int i = 0; i < len; i++) {
				result[i] /= norm;
			}
			return result;
		}

		public static int[][] FindChargePairs(int nsilac, SilacType silacType, double matchPpm,
		                                      double silacTimeCorrelationThreshold, SilacCluster[] silacClusters,
		                                      IsotopeCluster[] isotopeClusters, IPeakList peakList, float[] intensities,
		                                      double[] centerMz) {
			double[] m = new double[nsilac];
			for (int i = 0; i < m.Length; i++) {
				SilacCluster si = silacClusters[i];
				m[i] = GetUncalibratedSilacMass(si, silacType, isotopeClusters, centerMz, intensities);
			}
			int[] o = ArrayUtil.Order(m);
			m = ArrayUtil.SubArray(m, o);
			List<int[]> pairs = new List<int[]>();
			for (int i = 0; i < m.Length; i++) {
				SilacCluster si = silacClusters[o[i]];
				int chargei = isotopeClusters[si.IsotopeClusterindex1].Charge;
				double m1 = m[i];
				int ind1 = ArrayUtil.CeilIndex(m, m1 - matchPpm * 1e-6 * m1);
				for (int j = ind1; j < i; j++) {
					SilacCluster sj = silacClusters[o[j]];
					int chargej = isotopeClusters[sj.IsotopeClusterindex1].Charge;
					if (chargei == chargej) {
						continue;
					}
					int mini;
					double[] pi = CalcSilacClusterProfile(si, out mini, silacType, isotopeClusters, peakList,
					                                      intensities);
					int minj;
					double[] pj = CalcSilacClusterProfile(sj, out minj, silacType, isotopeClusters, peakList,
					                                      intensities);
					if (T03SilacAssembly.Correlate(pi, mini, pj, minj) > silacTimeCorrelationThreshold) {
						pairs.Add(new int[] {o[i], o[j]});
					}
				}
			}
			return pairs.ToArray();
		}

		private static double GetUncalibratedSilacMass(SilacCluster sc, SilacType type, IsotopeCluster[] isotopeClusters,
		                                               double[] centerMz, float[] intensities) {
			return
				GetUncalibratedFullIsotopePatternMassEstimate(sc.GetIsotopeClusterIndices(type), sc.GetIsotopePatternStarts(type),
				                                              sc.GetMassDiffs(type), isotopeClusters, centerMz, intensities);
		}

		private static double GetUncalibratedFullIsotopePatternMassEstimate(int[] isoClusterIndices, int[] isotopePatternStart,
		                                                                    double[] massShifts,
		                                                                    IsotopeCluster[] isotopeClusters,
		                                                                    double[] centerMz, float[] intensities) {
			double result = 0;
			double norm = 0;
			for (int i = 0; i < isoClusterIndices.Length; i++) {
				if (isoClusterIndices[i] < 0) {
					continue;
				}
				double intensity;
				double m =
					GetUncalibratedFullIsotopePatternMassEstimate(isotopeClusters[isoClusterIndices[i]], isotopePatternStart[i],
					                                              out intensity, centerMz, intensities);
				result += (m - massShifts[i]) * intensity;
				norm += intensity;
			}
			return result / norm;
		}

		private static double GetUncalibratedFullIsotopePatternMassEstimate(IsotopeCluster ic, int isotopePatternStart,
		                                                                    out double intensity, double[] centerMz,
		                                                                    float[] intensities) {
			double result = 0;
			intensity = 0;
			int charge = ic.Charge;
			for (int i = 0; i < ic.Count; i++) {
				if (i + isotopePatternStart >= 0) {
					int ind = ic.Members[i];
					double mz = centerMz[ind];
					double m = (mz - MolUtil.MassProton) * charge;
					result += intensities[ind] *
					          (m - MolUtil.GetAverageDifferenceToMonoisotope(m, i + isotopePatternStart));
					intensity += intensities[ind];
				}
			}
			if (intensity == 0) {
				for (int i = 0; i < ic.Count; i++) {
					int ind = ic.Members[i];
					double mz = centerMz[ind];
					double m = (mz - MolUtil.MassProton) * charge;
					result += intensities[ind] *
					          (m - (i + isotopePatternStart) * MolUtil.GetAverageDifferenceToMonoisotope(m, 1));
					intensity += intensities[ind];
				}
			}
			return result / intensity;
		}
	}
}