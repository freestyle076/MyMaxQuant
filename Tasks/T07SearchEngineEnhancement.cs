/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.IO;
using MaxQuant.Mol;
using MaxQuant.Num;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Tasks {
	public static class T07SearchEngineEnhancement {
		public static bool IsReverseProtein(string proteinId, string prefix) {
			return proteinId.IndexOf(prefix) >= 0;
		}

		public static Modification[] GetModifiedVersions(char aa, SilacType type, SilacLabel[] labels1, SilacLabel[] labels2) {
			AminoAcid a = AminoAcid.FromLetter(aa);
			switch (type) {
				case SilacType.Singlets:
					return new Modification[1];
				case SilacType.Doublets:
					Modification[] result = new Modification[2];
					if (a.HasFittingLabel(labels1)) {
						SilacLabel sl = a.GetFittingLabel(labels1);
						result[1] = Tables.modifications[AminoAcid.GetMascotModificationStringForLabel(sl)];
					}
					return result;
				case SilacType.Triplets:
					result = new Modification[3];
					if (a.HasFittingLabel(labels1)) {
						SilacLabel sl = a.GetFittingLabel(labels1);
						result[1] = Tables.modifications[AminoAcid.GetMascotModificationStringForLabel(sl)];
					}
					if (a.HasFittingLabel(labels2)) {
						SilacLabel sl = a.GetFittingLabel(labels2);
						result[2] = Tables.modifications[AminoAcid.GetMascotModificationStringForLabel(sl)];
					}
					return result;
			}
			throw new Exception("Never get here.");
		}

		public static bool ValidPeptide(string pep) {
			for (int i = 0; i < pep.Length; i++) {
				if (AminoAcid.singleLetterAas.IndexOf(pep[i]) == -1) {
					return false;
				}
			}
			return true;
		}

		public static AminoAcid[] CalcAllAas(SilacType type, SilacLabel[] labels1, SilacLabel[] labels2) {
			switch (type) {
				case SilacType.Singlets:
					return new AminoAcid[0];
				case SilacType.Doublets:
					return T08ReQuantification.CalcAllAas(labels1);
				case SilacType.Triplets:
					return T08ReQuantification.CalcAllAas(labels1, labels2);
			}
			throw new Exception("Never get here.");
		}

		public static double FindReverseHitThresholdValue(Identifications identifications, double totalPeptideFDR,
														  string reverseStr, IProteinSet proteinSet) {
			int n = identifications.Count;
			if (n == 0) {
				return 1;
			}
			double[] peps = new double[n];
			bool[] correct = new bool[n];
			for (int i = 0; i < n; i++) {
				MascotPeptide p = identifications.GetPeptidesAt(i)[0];
				peps[i] = p.Pep;
				correct[i] = !p.HasOnlyReverseHits(reverseStr, proteinSet);
			}
			int[] o = ArrayUtil.Order(peps);
			double forwardCount = 0;
			List<double> validPeps = new List<double>();
			for (int i = 0; i < n; i++) {
				int index = o[i];
				if (correct[index]) {
					forwardCount++;
				}
				double reverseCount = (i + 1) - forwardCount;
				if (reverseCount / forwardCount <= totalPeptideFDR) {
					validPeps.Add(peps[index]);
				}
			}
			if (validPeps.Count > 0) {
				return ArrayUtil.Max(validPeps.ToArray());
			}
			return 0;
		}

		public static void SetPep(Identifications identifications, BayesianInversion2D bi) {
			int n = identifications.Count;
			for (int i = 0; i < n; i++) {
				MascotPeptide p = identifications.GetPeptidesAt(i)[0];
				p.Pep = bi.GetValue(p.AltScore, Math.Log(p.Sequence.Length));
			}
			identifications.Write();
		}

		public static void FilterBySequence(Identifications identifications, HashSet<string> sequences) {
			int n = identifications.Count;
			List<int> oo = new List<int>();
			for (int i = 0; i < n; i++) {
				MascotPeptide p = identifications.GetPeptidesAt(i)[0];
				if (sequences.Contains(p.Sequence)) {
					oo.Add(i);
				}
			}
			identifications.Extract(oo.ToArray());
			identifications.Write();
		}

		public static void CollectSequences(Identifications identifications, double peptidePep,
		                                    HashSet<string> sequences) {
			int n = identifications.Count;
			for (int i = 0; i < n; i++) {
				MascotPeptide p = identifications.GetPeptidesAt(i)[0];
				if ((p.Pep <= peptidePep) || (double.IsNaN(p.Pep) && double.IsNaN(peptidePep))) {
					sequences.Add(p.Sequence);
				}
			}
		}

		public static void ApplyReverseHitThreshold(Identifications identifications, double pepThreshold) {
			int n = identifications.Count;
			List<int> o = new List<int>();
			for (int i = 0; i < n; i++) {
				MascotPeptide p = identifications.GetPeptidesAt(i)[0];
				if (p.Pep <= pepThreshold) {
					o.Add(i);
				}
			}
			identifications.Extract(o.ToArray());
			identifications.Write();
		}

		public static void FdrThresholding(string[] rawFiles, string decoyPrefix, double peptideFdr, double peptidePep,
								   MascotQueryType type, bool perFileThreshold, bool debug, HashSet<string> sequences,
								   bool onlyCollectSequences, IProteinSet proteinSet,
								   HashSet<string> labelModificationSet, IIdentificationProvider ip) {
			CalcFdr(rawFiles, decoyPrefix, type, proteinSet, debug, labelModificationSet, ip);
			double pepThreshVal;
			double[] pepThreshVals = new double[rawFiles.Length];
			if (peptideFdr < 1) {
				for (int i = 0; i < rawFiles.Length; i++) {
					pepThreshVals[i] = FindReverseHitThresholdValue(ip.GetIdentifications(rawFiles[i], type), peptideFdr,
																							   decoyPrefix, proteinSet);
					ip.Dispose();
					if (peptidePep<1){
						pepThreshVals[i] = Math.Min(pepThreshVals[i], peptidePep);
					}
				}
				pepThreshVal =ArrayUtil.Median(pepThreshVals);
				if (peptidePep < 1){
					pepThreshVal = Math.Min(pepThreshVal, peptidePep);
				}
			} else {
				for (int i = 0; i < rawFiles.Length; i++) {
					pepThreshVals[i] = peptidePep;
				}
				pepThreshVal = peptidePep;
			}
			if (onlyCollectSequences) {
				for (int i = 0; i < rawFiles.Length; i++) {
					CollectSequences(ip.GetIdentifications(rawFiles[i], type),
																perFileThreshold ? pepThreshVals[i] : pepThreshVal, sequences);
					ip.Dispose();
				}
			} else {
				for (int i = 0; i < rawFiles.Length; i++) {
					ApplyReverseHitThreshold(ip.GetIdentifications(rawFiles[i], type),
																		perFileThreshold ? pepThreshVals[i] : pepThreshVal);
					ip.Dispose();
				}
			}
		}

		public static void CalcFdr(string[] rawFiles, string decoyPrefix, MascotQueryType type, IProteinSet proteinSet, bool debug,
						HashSet<string> labelModificationSet, IIdentificationProvider ip) {
			CalcFdr(rawFiles, decoyPrefix, type, proteinSet, debug, ip);
		}

		private static void CalcFdr(string[] rawFiles, string decoyPrefix, MascotQueryType type, IProteinSet proteinSet, bool debug,
			IIdentificationProvider ip) {
			List<bool> correct = new List<bool>();
			List<double> scores = new List<double>();
			List<double> seqLen = new List<double>();
			for (int i = 0; i < rawFiles.Length; i++) {
				Identifications ident = ip.GetIdentifications(rawFiles[i], type);
				int n = ident.Count;
				for (int j = 0; j < n; j++) {
					bool c = ident.IsHighestScoringCorrect(j, decoyPrefix, proteinSet);
					double s = ident.GetHighestAltScore(j);
					double l = Math.Log(ident.GetBestSequence(j).Length);
					if (!double.IsNaN(s) && !double.IsInfinity(s)) {
						correct.Add(c);
						scores.Add(s);
						seqLen.Add(l);
					}
				}
				ip.Dispose();
				if (correct.Count > 10000000) {
					break;
				}
			}
			if (correct.Count == 0) {
				return;
			}
			bool write = debug && (type == MascotQueryType.Silac);
			BayesianInversion2D bi = new BayesianInversion2D(scores.ToArray(), seqLen.ToArray(), correct.ToArray(), write);
			if (write) {
				Write(rawFiles, bi);
			}
			for (int i = 0; i < rawFiles.Length; i++) {
				SetPep(ip.GetIdentifications(rawFiles[i], type), bi);
				ip.Dispose();
			}
		}

		private static void Write(string[] rawFiles, BayesianInversion2D bi) {
			string combinedFolder = rawFiles[0].Substring(0, rawFiles[0].LastIndexOf("\\")) + "\\combined";
			for (int len = 6; len <= 80; len++) {
				double loglen = Math.Log(len);
				string filename = "scoreDist" + len + ".txt";
				StreamWriter writer = new StreamWriter(combinedFolder + "\\" + filename);
				for (double score = 0; score < 600; score += 0.5) {
					double pep = bi.GetValue(score, loglen);
					double forw = bi.GetForwardHist(score, loglen);
					double reve = bi.GetReverseHist(score, loglen);
					writer.WriteLine(score + "\t" + pep + "\t" + forw + "\t" + reve);
				}
				writer.Close();
			}
		}

		public static void FilterBySequence(string[] rawFiles, string decoyPrefix, double totalPeptideFPR,
									double specificPeptideFPR, MascotQueryType type, bool perFileThreshold, bool debug,
									HashSet<string> sequences, bool onlyCollectSequences, IIdentificationProvider ip) {
			for (int i = 0; i < rawFiles.Length; i++) {
				FilterBySequence(ip.GetIdentifications(rawFiles[i], type), sequences);
				ip.Dispose();
			}
		}

		public static void FdrThresholding(string[] rawFiles, string[] recalFiles, string[] nonRecalFiles, string revstring,
								   double peptideFdr, double peptidePep, bool perFileThreshold, bool keepLowScorers,
								   IProteinSet proteinSet, HashSet<string> labelModificationSet, IIdentificationProvider ip, bool writeOut) {
			HashSet<string> sequences = new HashSet<string>();
			if (recalFiles.Length > 0) {
				FdrThresholding(recalFiles, revstring, peptideFdr, peptidePep, MascotQueryType.Silac,
														   perFileThreshold, writeOut, sequences, keepLowScorers, proteinSet,
														   labelModificationSet, ip);
				FdrThresholding(recalFiles, revstring, peptideFdr, peptidePep, MascotQueryType.Isotope,
														   perFileThreshold, false, sequences, keepLowScorers, proteinSet,
														   labelModificationSet, ip);
			}
			if (nonRecalFiles.Length > 0) {
				FdrThresholding(nonRecalFiles, revstring, peptideFdr, peptidePep, MascotQueryType.Silac,
														   perFileThreshold, false, sequences, keepLowScorers, proteinSet,
														   labelModificationSet, ip);
				FdrThresholding(nonRecalFiles, revstring, peptideFdr, peptidePep, MascotQueryType.Isotope,
														   perFileThreshold, false, sequences, keepLowScorers, proteinSet,
														   labelModificationSet, ip);
			}
			FdrThresholding(rawFiles, revstring, peptideFdr, peptidePep, MascotQueryType.Peak,
													   perFileThreshold,
													   false, sequences, keepLowScorers, proteinSet,
													   labelModificationSet, ip);
			if (keepLowScorers) {
				FilterBySequence(rawFiles, revstring, peptideFdr, peptidePep, MascotQueryType.Silac,
															perFileThreshold, writeOut, sequences, keepLowScorers, ip);
				FilterBySequence(rawFiles, revstring, peptideFdr, peptidePep, MascotQueryType.Isotope,
															perFileThreshold, false, sequences, keepLowScorers, ip);
				FilterBySequence(rawFiles, revstring, peptideFdr, peptidePep, MascotQueryType.Peak,
															perFileThreshold,
															false, sequences, keepLowScorers, ip);
			}
			LimitPep(rawFiles, MascotQueryType.Silac, ip);
			LimitPep(rawFiles, MascotQueryType.Isotope, ip);
			LimitPep(rawFiles, MascotQueryType.Peak, ip);
		}

		private static void LimitPep(string[] rawFiles, MascotQueryType type, IIdentificationProvider ip) {
			for (int i = 0; i < rawFiles.Length; i++) {
				Identifications ident = ip.GetIdentifications(rawFiles[i], type);
				int n = ident.Count;
				for (int j = 0; j < n; j++) {
					MascotPeptide p = ident.GetPeptidesAt(j)[0];
					p.Pep = Math.Min(p.Pep, 1);
				}
				ident.Write();
				ip.Dispose();
			}
		}
	}
}
