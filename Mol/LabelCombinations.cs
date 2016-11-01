/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using MaxQuant.Num;
using MaxQuant.Util;

namespace MaxQuant.Mol {
	public class LabelCombinations {
		private readonly double[][][] isotopeDistr1;
		private readonly double[][][] isotopeDistr2;
		private readonly Molecule[] labelingDiff1;
		private readonly Molecule[] labelingDiff2;
		private readonly double maxMassDiff;
		private readonly int[][] partitions;

		public LabelCombinations(SilacLabel[] labels, int maxNum, char[] allAaLetts)
			: this(new SilacLabel[0], labels, maxNum, allAaLetts) {
		}

		public LabelCombinations(SilacLabel[] labelsLow, SilacLabel[] labelsHigh, int maxNum, char[] allAaLetts) {
			AminoAcid[] aasHigh = new AminoAcid[labelsHigh.Length];
			char[] aaLettsHigh = new char[labelsHigh.Length];
			Molecule[] labelingDiff1High = new Molecule[labelsHigh.Length];
			Molecule[] labelingDiff2High = new Molecule[labelsHigh.Length];
			for (int i = 0; i < labelsHigh.Length; i++) {
				aasHigh[i] = AminoAcid.GetAminoAcidFromLabel(labelsHigh[i]);
				aaLettsHigh[i] = aasHigh[i].Letter;
				labelingDiff1High[i] = aasHigh[i].GetLabelingDiff1(labelsHigh[i]);
				labelingDiff2High[i] = aasHigh[i].GetLabelingDiff2(labelsHigh[i]);
			}
			AminoAcid[] aasLow = new AminoAcid[labelsLow.Length];
			char[] aaLettsLow = new char[labelsLow.Length];
			Molecule[] labelingDiff1Low = new Molecule[labelsLow.Length];
			Molecule[] labelingDiff2Low = new Molecule[labelsLow.Length];
			for (int i = 0; i < labelsLow.Length; i++) {
				aasLow[i] = AminoAcid.GetAminoAcidFromLabel(labelsLow[i]);
				aaLettsLow[i] = aasLow[i].Letter;
				labelingDiff1Low[i] = aasLow[i].GetLabelingDiff1(labelsLow[i]);
				labelingDiff2Low[i] = aasLow[i].GetLabelingDiff2(labelsLow[i]);
			}
			int[] indLow = new int[allAaLetts.Length];
			int[] indHigh = new int[allAaLetts.Length];
			for (int i = 0; i < allAaLetts.Length; i++) {
				indLow[i] = ArrayUtil.IndexOf(aaLettsLow, allAaLetts[i]);
				indHigh[i] = ArrayUtil.IndexOf(aaLettsHigh, allAaLetts[i]);
			}
			labelingDiff1 = new Molecule[allAaLetts.Length];
			labelingDiff2 = new Molecule[allAaLetts.Length];
			for (int i = 0; i < allAaLetts.Length; i++) {
				if (indLow[i] >= 0 && indHigh[i] >= 0) {
					labelingDiff1[i] = Molecule.Sum(labelingDiff2Low[indLow[i]], labelingDiff1High[indHigh[i]]);
					labelingDiff2[i] = Molecule.Sum(labelingDiff1Low[indLow[i]], labelingDiff2High[indHigh[i]]);
					Molecule[] w = Molecule.GetDifferences(labelingDiff1[i], labelingDiff2[i]);
					labelingDiff1[i] = w[0];
					labelingDiff2[i] = w[1];
				} else if (indLow[i] >= 0 && indHigh[i] < 0) {
					labelingDiff1[i] = labelingDiff2Low[indLow[i]];
					labelingDiff2[i] = labelingDiff1Low[indLow[i]];
				} else if (indLow[i] < 0 && indHigh[i] >= 0) {
					labelingDiff1[i] = labelingDiff1High[indHigh[i]];
					labelingDiff2[i] = labelingDiff2High[indHigh[i]];
				} else {
					throw new Exception("Should not happen.");
				}
			}
			partitions = GetPartitions(maxNum, allAaLetts.Length);
			isotopeDistr1 = new double[partitions.Length][][];
			isotopeDistr2 = new double[partitions.Length][][];
			for (int i = 0; i < partitions.Length; i++) {
				isotopeDistr1[i] = CalcDiff1(partitions[i]).GetIsotopeDistribution(0.2);
				isotopeDistr2[i] = CalcDiff2(partitions[i]).GetIsotopeDistribution(0.2);
			}
			maxMassDiff = 0;
			for (int i = 0; i < allAaLetts.Length; i++) {
				double diff = Math.Abs(labelingDiff1[i].GetMostLikelyMass(0.2) -
				                       labelingDiff2[i].GetMostLikelyMass(0.2));
				if (diff > maxMassDiff) {
					maxMassDiff = diff;
				}
			}
			maxMassDiff *= maxNum;
		}

		public double MaxMassDiff {
			get { return maxMassDiff; }
		}

		public int PartitionCount {
			get { return partitions.Length; }
		}

		public Molecule CalcDiff1(int[] counts) {
			return Molecule.Sum(labelingDiff1, counts);
		}

		public Molecule CalcDiff2(int[] counts) {
			return Molecule.Sum(labelingDiff2, counts);
		}

		public double[][] GetIsotopeDistribution1(int index) {
			return isotopeDistr1[index];
		}

		public double[][] GetIsotopeDistribution2(int index) {
			return isotopeDistr2[index];
		}

		private static int[][] GetPartitions(int maxNum, int nLabels) {
			int[][] p = NumUtil.GetPartitions(maxNum, nLabels + 1);
			for (int i = 0; i < p.Length; i++) {
				p[i] = ArrayUtil.SubArray(p[i], nLabels);
			}
			List<int[]> result = new List<int[]>();
			for (int i = 0; i < p.Length; i++) {
				if (ArrayUtil.Sum(p[i]) > 0) {
					result.Add(p[i]);
				}
			}
			return result.ToArray();
		}

		public int[] GetPartition(int i) {
			return partitions[i];
		}
	}
}