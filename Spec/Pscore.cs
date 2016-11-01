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

namespace MaxQuant.Spec {
	public static class Pscore {
		public static double CalcPtmScore(double ms2Tol, string ms2TolUnit, int topx, string sequence, Modification[] fixedModifications,
		                                  Modification[] variableModifications, int[] modCount, double[] specMasses, float[] specIntensities, double mz,
		                                  int charge, out bool incomplete, out PeptideModificationState outMods, out int counts,
		                                  out double delta, out MsmsPeakAnnotation[] description, out float[] intensities,
		                                  out float[] massDiffs, out Dictionary<ushort, double[]> modProb,
		                                  out Dictionary<ushort, double[]> modScoreDiffs, bool alternative) {
			if (!AminoAcid.ValidSequence(sequence)) {
				incomplete = true;
				outMods = null;
				counts = 0;
				delta = 0;
				description = null;
				intensities = null;
				massDiffs = null;
				modProb = null;
				modScoreDiffs = null;
				return 0;
			}
			Peptide p = new Peptide(sequence);
			p.ApplyFixedModifications(fixedModifications);
			PeptideModificationState[] mpeps = p.ApplyVariableModificationsFixedNumbers(variableModifications, modCount, out incomplete);
			counts = mpeps.Length;
			double maxScore = -double.MaxValue;
			double maxScore2 = -double.MaxValue;
			PeptideModificationState bestPep = null;
			description = null;
			intensities = null;
			massDiffs = null;
			double sumScore = 0;
			modProb = new Dictionary<ushort, double[]>();
			Dictionary<ushort, double[,]> scoreDiffs = new Dictionary<ushort, double[,]>();
			ushort[] w = mpeps[0].GetModifications();
			HashSet<ushort> um = new HashSet<ushort>();
			foreach (ushort x in w) {
				if (x != ushort.MaxValue && !um.Contains(x)) {
					um.Add(x);
				}
			}
			ushort[] uniqueMods = um.ToArray();
			foreach (ushort m in uniqueMods) {
				if (!modProb.ContainsKey(m)) {
					modProb.Add(m, new double[w.Length]);
					scoreDiffs.Add(m, new double[2,w.Length]);
				}
			}
			foreach (PeptideModificationState mpep in mpeps) {
				MsmsPeakAnnotation[] desc;
				float[] intens;
				float[] dm;
				double score = Score(specMasses, specIntensities, p, ms2Tol, ms2TolUnit, mz, charge,
				                     out desc, out intens, out dm, topx, mpep);
				ushort[] mm = mpep.GetModifications();
				double x = Math.Exp((score - 100) * NumUtil.log10 / 10.0);
				for (int i = 0; i < mm.Length; i++) {
					if (mm[i] != ushort.MaxValue) {
						modProb[mm[i]][i] += x;
					}
					foreach (ushort mod in uniqueMods) {
						if (mm[i] == mod) {
							if (score > scoreDiffs[mod][0, i]) {
								scoreDiffs[mod][0, i] = score;
							}
						} else {
							if (score > scoreDiffs[mod][1, i]) {
								scoreDiffs[mod][1, i] = score;
							}
						}
					}
				}
				sumScore += x;
				if (score > maxScore) {
					maxScore2 = maxScore;
					maxScore = score;
					bestPep = mpep;
					description = desc;
					intensities = intens;
					massDiffs = dm;
				} else if (score > maxScore2) {
					maxScore2 = score;
				}
			}
			modScoreDiffs = new Dictionary<ushort, double[]>();
			foreach (ushort mod in uniqueMods) {
				Modification modi = Tables.modificationList[mod];
				modScoreDiffs.Add(mod, new double[w.Length]);
				for (int i = 0; i < sequence.Length; i++) {
					char aa = sequence[i];
					if (modi.HasAA(aa)) {
						modScoreDiffs[mod][i] = scoreDiffs[mod][0, i] - scoreDiffs[mod][1, i];
					} else {
						modScoreDiffs[mod][i] = double.NaN;
					}
				}
			}
			foreach (double[] y in modProb.Values) {
				for (int i = 0; i < y.Length; i++) {
					y[i] /= sumScore;
				}
			}
			if (maxScore2 == -double.MaxValue) {
				delta = maxScore;
			} else {
				delta = maxScore - maxScore2;
			}
			outMods = bestPep;
			return maxScore;
		}

		public static double Score(double[] masses, double[] specMasses, float[] specIntensities, double tol, string tolUnit, double lnp, double lnq) {
			MsmsPeakAnnotation[] dummy;
			float[] intensities;
			float[] massDiffs;
			return Score(masses, null, specMasses, specIntensities, tol, tolUnit, out dummy, out intensities, out massDiffs, lnp, lnq);
		}

		public static double Score(double[] masses, MsmsPeakAnnotation[] descriptions, double[] specMasses, float[] specIntensities,
		                           double tol, string tolUnit, out MsmsPeakAnnotation[] description, out float[] intensities,
		                           out float[] massDiffs,
		                           double lnp, double lnq) {
			int n;
			int k =
				CountMatches(masses, descriptions, specMasses, specIntensities, tol, tolUnit, out description, out intensities, out massDiffs, out n);
			n = Math.Max(n, k);
			double s = -k * lnp - (n - k) * lnq -
			           NumericalRecipes.Factln(n) + NumericalRecipes.Factln(k) + NumericalRecipes.Factln(n - k);
			return 10.0 * s / NumUtil.log10;
		}

		public static double[] GetCidQueryMasses(string sequence, PeptideModificationState fixedMods,
		                                         PeptideModificationState varMods, double mz, int charge,
		                                         out MsmsPeakAnnotation[] description, bool includeNeutralLoss) {
			List<double> masses = new List<double>();
			List<MsmsPeakAnnotation> descriptions = new List<MsmsPeakAnnotation>();
			double[] y = GetYSeries(sequence, fixedMods, varMods);
			for (int i = 0; i < y.Length; i++) {
				double yval = y[i];
				masses.Add(yval);
				descriptions.Add(new MsmsPeakAnnotation(MsmsSeriesType.Y, i + 1, 1, false));
			}
			if (includeNeutralLoss) {
				double[] yn = GetNeutralLossYSeries(sequence, fixedMods, varMods);
				for (int i = 1; i < yn.Length; i++) {
					if (yn[i] > 0) {
						masses.Add(yn[i]);
						descriptions.Add(new MsmsPeakAnnotation(MsmsSeriesType.Y, i + 1, 1, true));
					}
				}
			}
			if (charge > 2) {
				for (int i = 0; i < y.Length; i++) {
					if (y[i] > 700) {
						double yval = y[i];
						double y2 = (yval + MolUtil.MassProton) / 2;
						masses.Add(y2);
						descriptions.Add(new MsmsPeakAnnotation(MsmsSeriesType.Y, i + 1, 2, false));
					}
				}
			}
			MsmsPeakAnnotation[] dx;
			double[] x = GetPhosphoXSeries(sequence, fixedMods, varMods, out dx, includeNeutralLoss);
			for (int i = 0; i < x.Length; i++) {
				masses.Add(x[i]);
				descriptions.Add(dx[i]);
			}
			double[] a = GetASeries(sequence, fixedMods, varMods);
			if (a.Length > 1) {
				masses.Add(a[1]);
				descriptions.Add(new MsmsPeakAnnotation(MsmsSeriesType.A, 2, 1, false));
			}
			double[] b = GetBSeries(sequence, fixedMods, varMods);
			for (int i = 1; i < b.Length; i++) {
				masses.Add(b[i]);
				descriptions.Add(new MsmsPeakAnnotation(MsmsSeriesType.B, i + 1, 1, false));
			}
			if (includeNeutralLoss) {
				double[] bn = GetNeutralLossBSeries(sequence, fixedMods, varMods);
				for (int i = 1; i < bn.Length; i++) {
					if (bn[i] > 0) {
						masses.Add(bn[i]);
						descriptions.Add(new MsmsPeakAnnotation(MsmsSeriesType.B, i + 1, 1, true));
					}
				}
			}
			description = descriptions.ToArray();
			return masses.ToArray();
		}

		public static double[] GetNeutralLossXSeries(string sequence, PeptideModificationState fixedMods,
		                                             PeptideModificationState varMods) {
			double[] result = new double[sequence.Length - 1];
			int neutralLossIndex = -1;
			for (int i = result.Length - 1; i >= 0; i--) {
				if (varMods.GetModificationAt(i) != ushort.MaxValue) {
					if (Tables.modificationList[varMods.GetModificationAt(i)].HasNeutralLoss) {
						neutralLossIndex = i;
						break;
					}
				}
			}
			if (neutralLossIndex != -1) {
				double[] x = GetXSeries(sequence, fixedMods, varMods);
				Modification mod = Tables.modificationList[varMods.GetModificationAt(neutralLossIndex)];
				double loss = mod.GetNeutralLoss()[0];
				for (int i = neutralLossIndex; i >= 0; i--) {
					result[i] = x[i] - loss;
				}
			}
			return result;
		}

		public static double[] GetNeutralLossYSeries(string sequence, PeptideModificationState fixedMods,
		                                             PeptideModificationState varMods) {
			double[] result = new double[sequence.Length - 1];
			int neutralLossIndex = -1;
			for (int i = sequence.Length - 1; i >= 1; i--) {
				if (varMods.GetModificationAt(i) != ushort.MaxValue) {
					if (Tables.modificationList[varMods.GetModificationAt(i)].HasNeutralLoss) {
						if (Tables.modificationList[varMods.GetModificationAt(i)].IsPhosphorylation && sequence[i] == 'Y') {
							continue;
						}
						neutralLossIndex = sequence.Length - 1 - i;
						break;
					}
				}
			}
			if (neutralLossIndex != -1) {
				double[] y = GetYSeries(sequence, fixedMods, varMods);
				Modification mod = Tables.modificationList[varMods.GetModificationAt(sequence.Length - 1 - neutralLossIndex)];
				double loss = mod.GetNeutralLoss()[0];
				for (int i = neutralLossIndex; i < sequence.Length - 1; i++) {
					result[i] = y[i] - loss;
				}
			}
			return result;
		}

		public static double[] GetNeutralLossBSeries(string sequence, PeptideModificationState fixedMods,
		                                             PeptideModificationState varMods) {
			double[] result = new double[sequence.Length - 1];
			int neutralLossIndex = -1;
			for (int i = 0; i < sequence.Length - 1; i++) {
				if (varMods.GetModificationAt(i) != ushort.MaxValue) {
					if (Tables.modificationList[varMods.GetModificationAt(i)].HasNeutralLoss) {
						if (Tables.modificationList[varMods.GetModificationAt(i)].IsPhosphorylation && sequence[i] == 'Y') {
							continue;
						}
						neutralLossIndex = i;
						break;
					}
				}
			}

			if (neutralLossIndex != -1) {
				double[] b = GetBSeries(sequence, fixedMods, varMods);
				Modification mod = Tables.modificationList[varMods.GetModificationAt(neutralLossIndex)];
				double loss = mod.GetNeutralLoss()[0];
				for (int i = neutralLossIndex; i < sequence.Length - 1; i++) {
					result[i] = b[i] - loss;
				}
			}
			return result;
		}

		public static double[] GetASeries(string sequence, PeptideModificationState fixedMods,
		                                  PeptideModificationState varMods) {
			double[] masses = GetBSeries(sequence, fixedMods, varMods);
			for (int i = 0; i < masses.Length; i++) {
				masses[i] += MolUtil.AIonMassOffset - MolUtil.BIonMassOffset;
			}
			return masses;
		}

		public static double[] GetBSeries(string sequence, PeptideModificationState fixedMods,
		                                  PeptideModificationState varMods) {
			double m = MolUtil.BIonMassOffset;
			if (fixedMods.GetNTermModification() != ushort.MaxValue) {
				Modification mod = Tables.modificationList[fixedMods.GetNTermModification()];
				m += mod.DeltaMass;
			}
			if (varMods.GetNTermModification() != ushort.MaxValue) {
				Modification mod = Tables.modificationList[varMods.GetNTermModification()];
				m += mod.DeltaMass;
			}
			double[] masses = new double[sequence.Length - 1];
			for (int i = 0; i < masses.Length; i++) {
				m += AminoAcid.aaMonoMasses[sequence[i]];
				if (fixedMods.GetModificationAt(i) != ushort.MaxValue) {
					Modification mod = Tables.modificationList[fixedMods.GetModificationAt(i)];
					m += mod.DeltaMass;
				}
				if (varMods.GetModificationAt(i) != ushort.MaxValue) {
					Modification mod = Tables.modificationList[varMods.GetModificationAt(i)];
					m += mod.DeltaMass;
				}
				masses[i] = m;
			}
			return masses;
		}

		public static double[] GetPhosphoXSeries(string sequence, PeptideModificationState fixedMods,
		                                         PeptideModificationState varMods,
		                                         out MsmsPeakAnnotation[] descriptions, bool includeNeutralLoss) {
			List<double> masses = new List<double>();
			List<MsmsPeakAnnotation> descr = new List<MsmsPeakAnnotation>();
			double[] y = GetYSeries(sequence, fixedMods, varMods);
			for (int i = 0; i < y.Length; i++) {
				int index = sequence.Length - 1 - i;
				if (varMods.GetModificationAt(index) != ushort.MaxValue) {
					if (Tables.modificationList[varMods.GetModificationAt(index)].IsPhosphorylation) {
						double x = y[i] + MolUtil.MassO + MolUtil.MassC;
						masses.Add(x);
						descr.Add(new MsmsPeakAnnotation(MsmsSeriesType.X, i + 1, 1, false));
						if (sequence[index] != 'Y' && includeNeutralLoss) {
							double xn = x - 97.976896;
							masses.Add(xn);
							descr.Add(new MsmsPeakAnnotation(MsmsSeriesType.X, i + 1, 1, true));
						}
					}
				}
			}
			descriptions = descr.ToArray();
			return masses.ToArray();
		}

		public static double[] GetXSeries(string sequence, PeptideModificationState fixedMods,
		                                  PeptideModificationState varMods) {
			double[] masses = GetYSeries(sequence, fixedMods, varMods);
			for (int i = 0; i < masses.Length; i++) {
				masses[i] += MolUtil.MassO + MolUtil.MassC;
			}
			return masses;
		}


		public static double[] GetYSeries(string sequence, PeptideModificationState fixedMods,
		                                  PeptideModificationState varMods) {
			double m = MolUtil.YIonMassOffset;
			if (fixedMods.GetCTermModification() != ushort.MaxValue) {
				Modification mod = Tables.modificationList[fixedMods.GetCTermModification()];
				m += mod.DeltaMass;
			}
			if (varMods.GetCTermModification() != ushort.MaxValue) {
				Modification mod = Tables.modificationList[varMods.GetCTermModification()];
				m += mod.DeltaMass;
			}
			int n = sequence.Length - 1;
			double[] masses = new double[n];
			for (int i = 0; i < masses.Length; i++) {
				m += AminoAcid.aaMonoMasses[sequence[n - i]];
				if (fixedMods.GetModificationAt(n - i) != ushort.MaxValue) {
					Modification mod = Tables.modificationList[fixedMods.GetModificationAt(n - i)];
					m += mod.DeltaMass;
				}
				if (varMods.GetModificationAt(n - i) != ushort.MaxValue) {
					Modification mod = Tables.modificationList[varMods.GetModificationAt(n - i)];
					m += mod.DeltaMass;
				}
				masses[i] = m;
			}
			return masses;
		}

		public static int CountMatches(double[] massesArray, MsmsPeakAnnotation[] descriptionsArray, double[] specMasses, float[] specIntensities,
		                               double tol, string tolUnit, out MsmsPeakAnnotation[] description,
		                               out float[] intensities, out float[] massDiffs, out int n) {
			if (specMasses == null) {
				n = 0;
				description = new MsmsPeakAnnotation[0];
				intensities = new float[0];
				massDiffs = new float[0];
				return 0;
			}
			int npeaks = specMasses.Length;
			if (npeaks == 0) {
				n = 0;
				description = new MsmsPeakAnnotation[0];
				intensities = new float[0];
				massDiffs = new float[0];
				return 0;
			}
			int k = 0;
			n = 0;
			List<MsmsPeakAnnotation> desc = null;
			List<float> intens = null;
			List<float> dm = null;
			if (descriptionsArray != null) {
				desc = new List<MsmsPeakAnnotation>();
				intens = new List<float>();
				dm = new List<float>();
			}
			for (int i = 0; i < massesArray.Length; i++) {
				double m = massesArray[i];
				n++;
				int ind = ArrayUtil.ClosestIndex(specMasses, m);
				if (ind != -1 && Match(m, specMasses[ind], tol, tolUnit)) {
					k++;
					if (descriptionsArray != null) {
						desc.Add(descriptionsArray[i]);
						intens.Add(specIntensities[ind]);
						dm.Add((float)(m - specMasses[ind]));
					}
				}
			}
			description = null;
			intensities = null;
			massDiffs = null;
			if (descriptionsArray != null) {
				description = desc.ToArray();
				intensities = intens.ToArray();
				massDiffs = dm.ToArray();
			}
			return k;
		}

		public static double Score(double[] specMasses, float[] specIntensities, Peptide pep, double tol, string tolUnit, double mz,
		                           int charge, out MsmsPeakAnnotation[] annot, out float[] intensities, out float[] massDiffs,
		                           int topx, PeptideModificationState modifications) {
			MsmsPeakAnnotation[] dummy;
			double p = Math.Min(6.0 / 100.0, 0.5);
			double lnp = Math.Log(p);
			double lnq = Math.Log(1.0 - p);
			double[] masses = GetCidQueryMasses(pep.Sequence, pep.FixedModifications, modifications, mz, charge, out dummy,
			                                    false);


			MsmsPeakAnnotation[] annot1;
			float[] intens1;
			float[] dm1;
			double score1 = Score(masses, dummy, specMasses, specIntensities, tol, tolUnit, out annot1, out intens1, out dm1, lnp, lnq);
			masses = GetCidQueryMasses(pep.Sequence, pep.FixedModifications, modifications, mz, charge, out dummy,
			                           true);

			MsmsPeakAnnotation[] annot2;
			float[] intens2;
			float[] dm2;
			double score2 = Score(masses, dummy, specMasses, specIntensities, tol, tolUnit, out annot2, out intens2, out dm2, lnp, lnq);
			if (score2 > score1) {
				annot = annot2;
				intensities = intens2;
				massDiffs = dm2;
				return score2;
			}
			annot = annot1;
			intensities = intens1;
			massDiffs = dm1;
			return score1;
		}

		public static bool Match(double m, double p, double tol, string tolUnit) {
			if (tolUnit.Equals("Da")) {
				return Math.Abs(m - p) <= tol;
			}
			if (tolUnit.Equals("ppm")) {
				return Math.Abs(m - p) / m * 1e6 <= tol;
			}
			throw new Exception("Unknown mass tolerance unit.");
		}
	}
}