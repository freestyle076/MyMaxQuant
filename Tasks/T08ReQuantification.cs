/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using MaxQuant.Mol;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Tasks {
	public static class T08ReQuantification {
		public static double ReQuantifyDoublets(string sequence, int silacState, double[] deltaMass1,
		                                        IsotopeCluster c, int icInd, IPeakList peakList, IRawFile rawFile,
		                                        out double[] outIntens, out int[] counts, SpectrumCache spectrumCache,
		                                        SilacLabel[] labels1, int maxSilacAa, bool subtractBackground,
		                                        int backgroundSubtractionQuantile) {
			AminoAcid[] allAas = CalcAllAas(labels1);
			char[] allAaLetts = AminoAcid.GetSingleLetters(allAas);
			LabelCombinations labelCombinations1 = new LabelCombinations(labels1, maxSilacAa, allAaLetts);
			AminoAcid[] a = new AminoAcid[labels1.Length];
			counts = new int[a.Length];
			double delta = 0;
			for (int j = 0; j < a.Length; j++) {
				a[j] = AminoAcid.GetAminoAcidFromLabel(labels1[j]);
				counts[j] = GetAaCount(a[j].Letter, sequence);
				delta += deltaMass1[j] * counts[j];
			}
			if (AllZero(counts)) {
				counts = null;
				outIntens = null;
				return double.NaN;
			}
			if (silacState == 1) {
				delta *= -1;
			}
			int[] members = c.Members;
			int charge = c.Charge;
			double[] masses = new double[members.Length];
			float[][] intens = new float[members.Length][];
			int[][] scanIndices = new int[members.Length][];
			for (int i = 0; i < members.Length; i++) {
				masses[i] = charge * (peakList.GetMz(members[i]) - MolUtil.MassProton) + delta;
				Peak peak = peakList.GetPeakDiscard(members[i]);
				scanIndices[i] = peak.GetScanIndices();
				intens[i] = IntegrateTranslatedPeak(peak, delta / charge, spectrumCache, rawFile,
				                                    subtractBackground, backgroundSubtractionQuantile);
				List<int> valids = new List<int>();
				for (int j = 0; j < intens[i].Length; j++) {
					if (intens[i][j] > 0) {
						valids.Add(j);
					}
				}
				int[] val = valids.ToArray();
				intens[i] = ArrayUtil.SubArray(intens[i], val);
				scanIndices[i] = ArrayUtil.SubArray(scanIndices[i], val);
			}
			List<int> valids1 = new List<int>();
			for (int i = 0; i < members.Length; i++) {
				if (intens[i].Length > 0) {
					valids1.Add(i);
				}
			}
			int[] val1 = valids1.ToArray();
			masses = ArrayUtil.SubArray(masses, val1);
			intens = ArrayUtil.SubArray(intens, val1);
			scanIndices = ArrayUtil.SubArray(scanIndices, val1);
			Molecule diff2 = labelCombinations1.CalcDiff1(counts);
			Molecule diff1 = labelCombinations1.CalcDiff2(counts);
			double[][] d1;
			double[][] d2;
			outIntens = new double[2];
			if (silacState == 0) {
				d1 = diff1.GetIsotopeDistribution(0.2);
				d2 = diff2.GetIsotopeDistribution(0.2);
				outIntens[0] = peakList.GetIsotopeClusterIntensity(icInd);
				outIntens[1] = GetClusterIntensity(intens);
			} else {
				d2 = diff1.GetIsotopeDistribution(0.2);
				d1 = diff2.GetIsotopeDistribution(0.2);
				outIntens[1] = peakList.GetIsotopeClusterIntensity(icInd);
				outIntens[0] = GetClusterIntensity(intens);
			}
			Dictionary<int, double>[] w = peakList.Vectorize(c, d1, masses, intens, scanIndices, d2);
			double r = T03SilacAssembly.FitRatio(w[0], w[1]);
			if (silacState == 1) {
				r = 1.0 / r;
			}
			return r;
		}

		public static double[] ReQuantifyTriplets(string sequence, int silacState, double[] deltaMass1,
		                                          double[] deltaMass2, IsotopeCluster c, int icInd, IPeakList peakList,
		                                          IRawFile rawFile, out double[] outIntens, out int[] counts,
		                                          SpectrumCache spectrumCache, SilacLabel[] labels1, SilacLabel[] labels2,
		                                          int maxSilacAa, bool subtractBackground, int backgroundSubtractionQuantile) {
			AminoAcid[] allAas = CalcAllAas(labels1, labels2);
			char[] allAaLetts = AminoAcid.GetSingleLetters(allAas);
			counts = new int[allAaLetts.Length];
			for (int j = 0; j < allAaLetts.Length; j++) {
				counts[j] = GetAaCount(allAaLetts[j], sequence);
			}
			if (AllZero(counts)) {
				counts = null;
				outIntens = null;
				return null;
			}
			double deltaLow = 0;
			for (int j = 0; j < labels1.Length; j++) {
				int cc = GetAaCount(AminoAcid.GetAminoAcidFromLabel(labels1[j]).Letter, sequence);
				deltaLow += deltaMass1[j] * cc;
			}
			double deltaHigh = 0;
			for (int j = 0; j < labels2.Length; j++) {
				int cc = GetAaCount(AminoAcid.GetAminoAcidFromLabel(labels2[j]).Letter, sequence);
				deltaHigh += deltaMass2[j] * cc;
			}
			double delta1;
			double delta2;
			CalcDeltasForTriplets(silacState, deltaLow, deltaHigh, out delta1, out delta2);
			int[] members = c.Members;
			int charge = c.Charge;
			double[] masses1 = new double[members.Length];
			float[][] intens1 = new float[members.Length][];
			int[][] scanIndices1 = new int[members.Length][];
			double[] masses2 = new double[members.Length];
			float[][] intens2 = new float[members.Length][];
			int[][] scanIndices2 = new int[members.Length][];
			for (int i = 0; i < members.Length; i++) {
				double m = charge * (peakList.GetMz(members[i]) - MolUtil.MassProton);
				masses1[i] = m + delta1;
				masses2[i] = m + delta2;
				Peak peak = peakList.GetPeakDiscard(members[i]);
				scanIndices1[i] = peak.GetScanIndices();
				scanIndices2[i] = peak.GetScanIndices();
				intens1[i] = IntegrateTranslatedPeak(peak, delta1 / charge, spectrumCache, rawFile,
				                                     subtractBackground, backgroundSubtractionQuantile);
				intens2[i] = IntegrateTranslatedPeak(peak, delta2 / charge, spectrumCache, rawFile,
				                                     subtractBackground, backgroundSubtractionQuantile);
				List<int> valids = new List<int>();
				for (int j = 0; j < intens1[i].Length; j++) {
					if (intens1[i][j] > 0) {
						valids.Add(j);
					}
				}
				int[] val = valids.ToArray();
				intens1[i] = ArrayUtil.SubArray(intens1[i], val);
				scanIndices1[i] = ArrayUtil.SubArray(scanIndices1[i], val);
				valids = new List<int>();
				for (int j = 0; j < intens2[i].Length; j++) {
					if (intens2[i][j] > 0) {
						valids.Add(j);
					}
				}
				val = valids.ToArray();
				intens2[i] = ArrayUtil.SubArray(intens2[i], val);
				scanIndices2[i] = ArrayUtil.SubArray(scanIndices2[i], val);
			}
			List<int> valids1 = new List<int>();
			for (int i = 0; i < members.Length; i++) {
				if (intens1[i].Length > 0) {
					valids1.Add(i);
				}
			}
			int[] val1 = valids1.ToArray();
			masses1 = ArrayUtil.SubArray(masses1, val1);
			intens1 = ArrayUtil.SubArray(intens1, val1);
			scanIndices1 = ArrayUtil.SubArray(scanIndices1, val1);
			List<int> valids2 = new List<int>();
			for (int i = 0; i < members.Length; i++) {
				if (intens2[i].Length > 0) {
					valids2.Add(i);
				}
			}
			int[] val2 = valids2.ToArray();
			masses2 = ArrayUtil.SubArray(masses2, val2);
			intens2 = ArrayUtil.SubArray(intens2, val2);
			scanIndices2 = ArrayUtil.SubArray(scanIndices2, val2);
			LabelCombinations labelCombinations01 = new LabelCombinations(labels1, maxSilacAa, allAaLetts);
			LabelCombinations labelCombinations02 = new LabelCombinations(labels2, maxSilacAa, allAaLetts);
			Molecule x1 = labelCombinations01.CalcDiff1(counts);
			Molecule x2 = labelCombinations01.CalcDiff2(counts);
			Molecule y1 = labelCombinations02.CalcDiff1(counts);
			Molecule y2 = labelCombinations02.CalcDiff2(counts);
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
			double[][] d1 = diff1.GetIsotopeDistribution(0.2);
			double[][] d2 = diff2.GetIsotopeDistribution(0.2);
			double[][] d3 = diff3.GetIsotopeDistribution(0.2);
			double ratio10;
			double ratio20;
			double ratio21;
			outIntens = new double[3];
			switch (silacState) {
				case 0: {
					Dictionary<int, double>[] w10 = peakList.Vectorize(c, d1, masses1, intens1, scanIndices1, d2);
					ratio10 = T03SilacAssembly.FitRatio(w10[0], w10[1]);
					Dictionary<int, double>[] w20 = peakList.Vectorize(c, d1, masses2, intens2, scanIndices2, d3);
					ratio20 = T03SilacAssembly.FitRatio(w20[0], w20[1]);
					Dictionary<int, double>[] w21 =
						peakList.Vectorize(masses1, intens1, scanIndices1, d2, masses2, intens2, scanIndices2, d3, c.Charge);
					ratio21 = T03SilacAssembly.FitRatio(w21[0], w21[1]);
					outIntens[0] = peakList.GetIsotopeClusterIntensity(icInd);
					outIntens[1] = GetClusterIntensity(intens1);
					outIntens[2] = GetClusterIntensity(intens2);
					break;
				}
				case 1: {
					Dictionary<int, double>[] w01 = peakList.Vectorize(c, d2, masses1, intens1, scanIndices1, d1);
					ratio10 = 1.0 / T03SilacAssembly.FitRatio(w01[0], w01[1]);
					Dictionary<int, double>[] w20 =
						peakList.Vectorize(masses1, intens1, scanIndices1, d1, masses2, intens2, scanIndices2, d3, c.Charge);
					ratio20 = T03SilacAssembly.FitRatio(w20[0], w20[1]);
					Dictionary<int, double>[] w21 = peakList.Vectorize(c, d2, masses2, intens2, scanIndices2, d3);
					ratio21 = T03SilacAssembly.FitRatio(w21[0], w21[1]);
					outIntens[1] = peakList.GetIsotopeClusterIntensity(icInd);
					outIntens[0] = GetClusterIntensity(intens1);
					outIntens[2] = GetClusterIntensity(intens2);
					break;
				}
				case 2: {
					Dictionary<int, double>[] w10 =
						peakList.Vectorize(masses1, intens1, scanIndices1, d1, masses2, intens2, scanIndices2, d2, c.Charge);
					ratio10 = T03SilacAssembly.FitRatio(w10[0], w10[1]);
					Dictionary<int, double>[] w02 = peakList.Vectorize(c, d3, masses1, intens1, scanIndices1, d1);
					ratio20 = 1.0 / T03SilacAssembly.FitRatio(w02[0], w02[1]);
					Dictionary<int, double>[] w12 = peakList.Vectorize(c, d3, masses2, intens2, scanIndices2, d2);
					ratio21 = 1.0 / T03SilacAssembly.FitRatio(w12[0], w12[1]);
					outIntens[2] = peakList.GetIsotopeClusterIntensity(icInd);
					outIntens[0] = GetClusterIntensity(intens1);
					outIntens[1] = GetClusterIntensity(intens2);
					break;
				}
				default:
					throw new Exception("Impossible.");
			}
			return new double[] {ratio10, ratio20, ratio21};
		}

		public static float[] IntegrateTranslatedPeak(Peak peak, double delta, SpectrumCache spectrumCache,
		                                              IRawFile rawFile, bool subtractBackground,
		                                              int backgroundSubtractionQuantile) {
			int[] scanIndices = peak.GetScanIndices();
			float[] intensities = peak.GetOriginalIntensities();
			int q = ArrayUtil.MaxInd(intensities);
			float maxInt = intensities[q];
			int maxScanIndex = scanIndices[q];
			float[] result = new float[scanIndices.Length];
			for (int i = 0; i < scanIndices.Length; i++) {
				int scanIndex = scanIndices[i];
				if (intensities[i] < maxInt * 0.01) {
					continue;
				}
				if (scanIndex > maxScanIndex + 20 || scanIndex < maxScanIndex - 20) {
					continue;
				}
				Spectrum spectrum;
				if (!spectrumCache.ContainsScanIndex(scanIndex)) {
					spectrum = rawFile.GetMS1Spectrum(scanIndex, subtractBackground, backgroundSubtractionQuantile);
					spectrumCache.Add(scanIndex, spectrum);
				} else {
					spectrum = spectrumCache[scanIndex];
				}
				double minMass = peak.GetMinMass(i) + delta;
				double maxMass = peak.GetMaxMass(i) + delta;
				int minInd = spectrum.GetCeilIndex(minMass);
				int maxInd = spectrum.GetFloorIndex(maxMass);
				if (minInd == -1 || maxInd == -1) {
					result[i] = 0;
				} else {
					float[] intensityProfile = new float[maxInd - minInd + 1];
					for (int ind = minInd; ind <= maxInd; ind++) {
						intensityProfile[ind - minInd] = spectrum.GetIntensity(ind);
					}
					result[i] = IntegrateProfile(intensityProfile, T01PeakDetection.maxIntensity);
				}
			}
			return result;
		}

		private static float IntegrateProfile(float[] profile, bool maxIntensity) {
			if (!HasMaximum(profile)) {
				return 0;
			}
			float peakIntensity = 0;
			for (int j = 0; j < profile.Length; j++) {
				float intensity = profile[j];
				if (maxIntensity) {
					if (intensity > peakIntensity) {
						peakIntensity = intensity;
					}
				} else {
					peakIntensity += intensity;
				}
			}
			return peakIntensity;
		}

		private static bool HasMaximum(float[] profile) {
			for (int i = 1; i < profile.Length - 1; i++) {
				if (profile[i - 1] < profile[i] && profile[i + 1] < profile[i]) {
					return true;
				}
			}
			return false;
		}

		public static void CalcDeltasForTriplets(int silacState, double deltaLow, double deltaHigh, out double delta1,
		                                         out double delta2) {
			switch (silacState) {
				case 0:
					delta1 = deltaLow;
					delta2 = deltaHigh;
					break;
				case 1:
					delta1 = -deltaLow;
					delta2 = deltaHigh - deltaLow;
					break;
				case 2:
					delta1 = -deltaHigh;
					delta2 = -deltaHigh + deltaLow;
					break;
				default:
					throw new Exception("Impossible.");
			}
		}

		public static bool AllZero(int[] counts) {
			for (int i = 0; i < counts.Length; i++) {
				if (counts[i] != 0) {
					return false;
				}
			}
			return true;
		}

		public static float GetClusterIntensity(IEnumerable<float[]> intens) {
			float max = 0;
			foreach (float[] x in intens) {
				foreach (float y in x) {
					if (y > max) {
						max = y;
					}
				}
			}
			return max;
		}

		public static int GetAaCount(char aa, string sequence) {
			int c = 0;
			for (int i = 0; i < sequence.Length; i++) {
				if (sequence[i] == aa) {
					c++;
				}
			}
			return c;
		}

		public static AminoAcid[] CalcAllAas(SilacLabel[] labels1) {
			List<AminoAcid> result = new List<AminoAcid>();
			for (int i = 0; i < labels1.Length; i++) {
				foreach (AminoAcid aa in AminoAcid.aminoAcids) {
					if (aa.Equals(AminoAcid.GetAminoAcidFromLabel(labels1[i]))) {
						result.Add(aa);
						break;
					}
				}
			}
			return result.ToArray();
		}

		public static AminoAcid[] CalcAllAas(SilacLabel[] labels1, SilacLabel[] labels2) {
			List<AminoAcid> result = new List<AminoAcid>();
			foreach (AminoAcid aa in AminoAcid.aminoAcids) {
				bool found = false;
				for (int i = 0; i < labels1.Length; i++) {
					if (aa.Equals(AminoAcid.GetAminoAcidFromLabel(labels1[i]))) {
						result.Add(aa);
						found = true;
						break;
					}
				}
				if (!found) {
					for (int i = 0; i < labels2.Length; i++) {
						if (aa.Equals(AminoAcid.GetAminoAcidFromLabel(labels2[i]))) {
							result.Add(aa);
							break;
						}
					}
				}
			}
			return result.ToArray();
		}
	}
}