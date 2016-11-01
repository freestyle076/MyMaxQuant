/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Collections.Generic;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Tasks {
	public static class T10ProteinFdr {
		public static string[] AllPeptideSequences(IIdentifiedPeptideList peptides) {
			string[] result = new string[peptides.Count];
			for (int i = 0; i < result.Length; i++) {
				result[i] = peptides.GetSequence(i);
			}
			return result;
		}

		public static int[] ExtractProteinGroups(TmpIdentifiedProteinGroup[] tmpProteinGroups, IIdentifiedPeptideList peptides,
										  double proteinFDR, PeptideScoring peptideScoring, int nsilac, int minPeptides,
										  int minUniquePeptides, string reversePrefix) {
			double[] peps = new double[tmpProteinGroups.Length];
			string[] sequences = AllPeptideSequences(peptides);
			for (int i = 0; i < tmpProteinGroups.Length; i++) {
				peps[i] = tmpProteinGroups[i].GetPep(peptides, sequences, minPeptides, minUniquePeptides, peptideScoring, nsilac);
			}
			int[] o = ArrayUtil.Order(peps);
			int revCount = 0;
			List<double> goodThresholdValues = new List<double>();
			for (int i = 0; i < tmpProteinGroups.Length; i++) {
				if (tmpProteinGroups[o[i]].HasReverseProteins(reversePrefix)) {
					revCount++;
				}
				double forwCount = i + 1.0 - revCount;
				if (revCount / forwCount <= proteinFDR && peps[o[i]] < 1.01) {
					goodThresholdValues.Add(peps[o[i]]);
				}
			}
			double thresholdValue = ArrayUtil.Max(goodThresholdValues.ToArray());
			revCount = 0;
			List<int> oo = new List<int>();
			for (int i = 0; i < tmpProteinGroups.Length; i++) {
				if (tmpProteinGroups[o[i]].HasReverseProteins(reversePrefix)) {
					revCount++;
				}
				double forwCount = i + 1.0 - revCount;
				if (revCount / forwCount <= proteinFDR || peps[o[i]] < thresholdValue) {
					oo.Add(o[i]);
				}
			}
			string[] dummy1 = new string[oo.Count];
			for (int i = 0; i < oo.Count; i++) {
				dummy1[i] = tmpProteinGroups[oo[i]].ProteinIds[0];
			}
			int[] ooo = ArrayUtil.Order(dummy1);
			o = ArrayUtil.SubArray(oo, ooo);
			return o;
		}
	}
}
