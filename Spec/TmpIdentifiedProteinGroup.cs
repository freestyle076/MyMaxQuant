/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using MaxQuant.Spec;
using MaxQuant.Tasks;

namespace MaxQuant.Spec {
	public class TmpIdentifiedProteinGroup {
		private readonly string[] proteinIds;
		private readonly string[] peptideSequences;

		public TmpIdentifiedProteinGroup(string[] proteinIds, string[] peptideSequences) {
			this.proteinIds = proteinIds;
			this.peptideSequences = peptideSequences;
		}

		public string[] ProteinIds {
			get { return proteinIds; }
		}

		public string[] PeptideSequences {
			get { return peptideSequences; }
		}

		public double GetPep(IIdentifiedPeptideList peptides, string[] sequences, int minPeptides, int minUniquePeptides,
		                     PeptideScoring peptideScoring, int nsilac) {
			if(peptideSequences.Length < minPeptides) {
				return 1.1;
			}
			double result = 1;
			int uniqueCount = 0;
			for (int i = 0; i < peptideSequences.Length; i++) {
				int index = Array.BinarySearch(sequences, peptideSequences[i]);
				if(peptides.IsUniqueGroup(index)) {
					uniqueCount++;
				}
				result *= peptides.GetPep(index, peptideScoring, nsilac);
			}
			if(uniqueCount < minUniquePeptides) {
				return 1.1;
			}
			return result;
		}

		public bool HasReverseProteins(string revstring) {
			for (int i = 0; i < proteinIds.Length; i++) {
				if (T07SearchEngineEnhancement.IsReverseProtein(proteinIds[i], revstring)) {
					return true;
				}
			}
			return false;
		}
	}
}