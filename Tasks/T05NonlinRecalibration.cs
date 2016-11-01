/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using MaxQuant.Mol;
using MaxQuant.Num;
using MaxQuant.Spec;

namespace MaxQuant.Tasks {
	public static class T05NonlinRecalibration {
		public const int nparamsIntensityFit = 1;
		public const int nparamsMzFit = 5;

		private static readonly double[] startValsMzFit =
			new double[] {-0.320846816779201, -0.0563056487654107, 0.00433475125727255, 0.00069497234303954, 0, 0, 0};

		public static void RecalibrateImpl(SilacType silacType, int nranges, int[] iranges, int[][] chargePairs,
		                                   SilacCluster[] silacClusters, ref double[,] mzCalibrationParams,
		                                   ref double[,] intensityCalibrationParams, double massErrorNormalization,
		                                   IsotopeCluster[] isotopeCluster, IPeakList peakList, float[] intensities) {
			double[] x = new double[chargePairs.Length];
			double[] y = new double[chargePairs.Length];
			double[] sig = new double[chargePairs.Length];
			for (int i = 0; i < chargePairs.Length; i++) {
				x[i] = i;
				int ind1 = chargePairs[i][0];
				int ind2 = chargePairs[i][1];
				SilacCluster sc1 = silacClusters[ind1];
				SilacCluster sc2 = silacClusters[ind2];
				sig[i] = Math.Sqrt(sc1.GetMassError(massErrorNormalization) * sc1.GetMassError(massErrorNormalization) +
				                   sc2.GetMassError(massErrorNormalization) * sc2.GetMassError(massErrorNormalization));
			}
			intensityCalibrationParams = new double[nranges,nparamsIntensityFit];
			double[,] icp = intensityCalibrationParams;
			Funcs2 f = delegate(double x1, double[] aa) {
				int ind = (int) Math.Round(x1);
				int ind1 = chargePairs[ind][0];
				int ind2 = chargePairs[ind][1];
				double[,] mzAa = new double[nranges,nparamsMzFit];
				int c = 0;
				for (int j = 0; j < iranges.Length; j++) {
					for (int i = 0; i < nparamsMzFit; i++) {
						mzAa[iranges[j], i] = aa[c++];
					}
				}
				double y1 = GetCalibratedSilacMassStat(silacClusters[ind1], mzAa, icp, silacType, isotopeCluster, peakList,
				                                       intensities);
				double y2 = GetCalibratedSilacMassStat(silacClusters[ind2], mzAa, icp, silacType, isotopeCluster, peakList,
				                                       intensities);
				double yy = y1 - y2;
				return yy;
			};
			double chisq;
			double[] a = new double[nparamsMzFit * iranges.Length];
			int county = 0;
			for (int j = 0; j < iranges.Length; j++) {
				for (int i = 0; i < nparamsMzFit; i++) {
					a[county++] = startValsMzFit[i];
				}
			}
			NumUtil.FitNonlin(x, y, sig, a, out chisq, f);
			mzCalibrationParams = new double[nranges,nparamsMzFit];
			int s = 0;
			for (int j = 0; j < iranges.Length; j++) {
				for (int i = 0; i < nparamsMzFit; i++) {
					mzCalibrationParams[iranges[j], i] = a[s++];
				}
			}
			double[,] mcp = mzCalibrationParams;
			f = delegate(double x1, double[] aa) {
				int ind = (int) Math.Round(x1);
				int ind1 = chargePairs[ind][0];
				int ind2 = chargePairs[ind][1];
				double[,] intensityAa = new double[nranges,nparamsIntensityFit];
				int c = 0;
				for (int j = 0; j < iranges.Length; j++) {
					for (int i = 0; i < nparamsIntensityFit; i++) {
						intensityAa[iranges[j], i] = aa[c++];
					}
				}
				double y1 = GetCalibratedSilacMassStat(silacClusters[ind1], mcp, intensityAa, silacType, isotopeCluster, peakList,
				                                       intensities);
				double y2 = GetCalibratedSilacMassStat(silacClusters[ind2], mcp, intensityAa, silacType, isotopeCluster, peakList,
				                                       intensities);
				double yy = y1 - y2;
				return yy;
			};
			a = new double[nparamsIntensityFit * iranges.Length];
			NumUtil.FitNonlin(x, y, sig, a, out chisq, f);
			s = 0;
			for (int j = 0; j < iranges.Length; j++) {
				for (int i = 0; i < nparamsIntensityFit; i++) {
					intensityCalibrationParams[iranges[j], i] = a[s++];
				}
			}
			icp = intensityCalibrationParams;
			f = delegate(double x1, double[] aa) {
				int ind = (int) Math.Round(x1);
				int ind1 = chargePairs[ind][0];
				int ind2 = chargePairs[ind][1];
				double[,] mzAa = new double[nranges,nparamsMzFit];
				int c = 0;
				for (int j = 0; j < iranges.Length; j++) {
					for (int i = 0; i < nparamsMzFit; i++) {
						mzAa[iranges[j], i] = aa[c++];
					}
				}
				double y1 = GetCalibratedSilacMassStat(silacClusters[ind1], mzAa, icp, silacType, isotopeCluster, peakList,
				                                       intensities);
				double y2 = GetCalibratedSilacMassStat(silacClusters[ind2], mzAa, icp, silacType, isotopeCluster, peakList,
				                                       intensities);
				double yy = y1 - y2;
				return yy;
			};
			a = new double[nparamsMzFit * iranges.Length];
			county = 0;
			for (int j = 0; j < iranges.Length; j++) {
				for (int i = 0; i < nparamsMzFit; i++) {
					a[county++] = mzCalibrationParams[iranges[j], i];
				}
			}
			NumUtil.FitNonlin(x, y, sig, a, out chisq, f);
			s = 0;
			for (int j = 0; j < iranges.Length; j++) {
				for (int i = 0; i < nparamsMzFit; i++) {
					mzCalibrationParams[iranges[j], i] = a[s++];
				}
			}
		}

		private static double GetCalibratedSilacMassStat(SilacCluster sc, double[,] mzCalibrationPar,
		                                                 double[,] intensityCalibrationPar,
		                                                 SilacType type, IsotopeCluster[] isotopeCluster, IPeakList peakList,
		                                                 float[] intensities) {
			int dummy;
			return
				GetFullIsotopePatternMassEstimateNoBootstrapStat(sc.GetIsotopeClusterIndices(type), sc.GetIsotopePatternStarts(type),
				                                                 sc.GetMassDiffs(type), mzCalibrationPar, intensityCalibrationPar,
				                                                 true,
				                                                 out dummy, isotopeCluster, peakList, intensities);
		}

		private static double GetFullIsotopePatternMassEstimateNoBootstrapStat(int[] isoClusterIndices,
		                                                                       int[] isotopePatternStart,
		                                                                       double[] massShifts, double[,] mzCalibrationPar,
		                                                                       double[,] intensityCalibrationPar, bool discard,
		                                                                       out int npoints,
		                                                                       IsotopeCluster[] isotopeClusters,
		                                                                       IPeakList peakList, float[] intensities) {
			npoints = 0;
			double result = 0;
			double norm = 0;
			for (int i = 0; i < isoClusterIndices.Length; i++) {
				if (isoClusterIndices[i] < 0) {
					continue;
				}
				double intensity;
				int np;
				double m =
					GetFullIsotopePatternMassEstimateNoBootstrapStat(isotopeClusters[isoClusterIndices[i]], isotopePatternStart[i],
					                                                 out intensity,
					                                                 mzCalibrationPar, intensityCalibrationPar, discard, out np,
					                                                 peakList, intensities);
				result += (m - massShifts[i]) * intensity;
				norm += intensity;
				npoints += np;
			}
			return result / norm;
		}

		private static double GetFullIsotopePatternMassEstimateNoBootstrapStat(IsotopeCluster ic, int isotopePatternStart,
		                                                                       out double intensity,
		                                                                       double[,] mzCalibrationPar,
		                                                                       double[,] intensityCalibrationPar, bool discard,
		                                                                       out int npoints, IPeakList peakList,
		                                                                       float[] intensities) {
			double result = 0;
			intensity = 0;
			npoints = 0;
			int charge = ic.Charge;
			for (int i = 0; i < ic.Count; i++) {
				if (i + isotopePatternStart >= 0) {
					int ind = ic.Members[i];
					int np;
					double mz = discard
					            	? peakList.GetMzDiscard(ind, mzCalibrationPar, intensityCalibrationPar, out np)
					            	:
					            		peakList.GetMz(ind, mzCalibrationPar, intensityCalibrationPar, out np);
					double m = (mz - MolUtil.MassProton) * charge;
					result += intensities[ind] *
					          (m - MolUtil.GetAverageDifferenceToMonoisotope(m, i + isotopePatternStart));
					intensity += intensities[ind];
					npoints += np;
				}
			}
			if (intensity == 0) {
				for (int i = 0; i < ic.Count; i++) {
					int ind = ic.Members[i];
					int np;
					double mz = discard
					            	? peakList.GetMzDiscard(ind, mzCalibrationPar, intensityCalibrationPar, out np)
					            	: peakList.GetMz(ind, mzCalibrationPar, intensityCalibrationPar, out np);
					double m = (mz - MolUtil.MassProton) * charge;
					result += intensities[ind] *
					          (m - (i + isotopePatternStart) * MolUtil.GetAverageDifferenceToMonoisotope(m, 1));
					intensity += intensities[ind];
					npoints += np;
				}
			}
			return result / intensity;
		}
	}
}