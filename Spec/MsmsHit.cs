/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.IO;
using MaxQuant.Mol;
using MaxQuant.Spec;
using MaxQuant.Tasks;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public class MsmsHit {
		private readonly int rawFileIndex;
		private readonly int scanNumber;
		private readonly int silacIndex = -1;
		private readonly MascotQueryType type;
		private readonly double pep;
		private readonly double score;
		private readonly double altScore;
		private readonly double deltaScoreAll;
		private readonly double deltaScoreDifferentPep;
		private PeptideModificationState labelModifications;
		private PeptideModificationState trueModifications;
		private int id;
		private int evidenceId;
		private readonly double elutionTime;
		private readonly double mz;
		private readonly double monoisotopicMz;
		private double ptmScore;
		private double deltaPtmScore;
		private int ptmScoreCounts;
		private bool ptmScoreFinished;
		private readonly double massErrorPpm;
		private bool changed;
		private MsmsPeakAnnotation[] description;
		private float[] intensities;
		private float[] massDiffs;
		private readonly int charge;
		private Dictionary<ushort, double[]> modProbabilities;
		private Dictionary<ushort, double[]> modScoreDiffs;


		public MsmsHit(MascotPeptide[] peptides, int scanNumber, int fileIndex, MascotQueryType type,
		               int silacIndex, double elutionTime, double mz, double monoisotopicMz, ushort[] fixedMods,
		               HashSet<string> labelModificationSet, SilacType silacType, SilacLabel[] labels1, SilacLabel[] labels2) {
			rawFileIndex = fileIndex;
			this.scanNumber = scanNumber;
			this.silacIndex = silacIndex;
			this.type = type;
			pep = peptides[0].Pep;
			score = peptides[0].MascotScore;
			altScore = peptides[0].AltScore;
			charge = (int) Math.Round(peptides[0].Mass / mz);
			massErrorPpm = peptides[0].MassErrorPpm;
			if (peptides.Length > 1) {
				deltaScoreAll = peptides[0].MascotScore - peptides[1].MascotScore;
			} else {
				deltaScoreAll = peptides[0].MascotScore;
			}
			int ind = -1;
			for (int i = 1; i < peptides.Length; i++) {
				if (!peptides[i].Sequence.Equals(peptides[0].Sequence)) {
					ind = i;
					break;
				}
			}
			if (ind != -1) {
				deltaScoreDifferentPep = peptides[0].MascotScore - peptides[ind].MascotScore;
			} else {
				deltaScoreDifferentPep = peptides[0].MascotScore;
			}
			labelModifications =
				peptides[0].Modifications.GetLabelModifications(fixedMods, peptides[0].Sequence, labelModificationSet);
			trueModifications = peptides[0].Modifications.GetTrueModifications(labelModificationSet);
			this.elutionTime = elutionTime;
			this.mz = mz;
			this.monoisotopicMz = monoisotopicMz;
			if (silacType == SilacType.Singlets) {
				this.silacIndex = 0;
			} else if (type != MascotQueryType.Silac) {
				AminoAcid[] allAas = T07SearchEngineEnhancement.CalcAllAas(silacType, labels1, labels2);
				char[] allAaLetts = AminoAcid.GetSingleLetters(allAas);
				for (int i = 0; i < allAaLetts.Length; i++) {
					this.silacIndex =
						CalcSilacIndex(allAaLetts[i], labelModifications, peptides[0].Sequence, silacType, labels1,
						               labels2);
					if (this.silacIndex != -1) {
						break;
					}
				}
			}
		}

		private static int CalcSilacIndex(char aa, PeptideModificationState labelMods, string sequence, SilacType type,
		                                  SilacLabel[] labels1, SilacLabel[] labels2) {
			for (int i = 0; i < sequence.Length; i++) {
				if (sequence[i] != aa) {
					continue;
				}
				ushort m = labelMods.GetModificationAt(i);
				if (type == SilacType.Doublets) {
					if (m != ushort.MaxValue) {
						return 1;
					}
					return 0;
				}
				if (type == SilacType.Triplets) {
					int ind1 = -1;
					for (int j = 0; j < labels1.Length; j++) {
						if (AminoAcid.GetAminoAcidFromLabel(labels1[j]).Letter == aa) {
							ind1 = j;
						}
					}
					int ind2 = -1;
					for (int j = 0; j < labels2.Length; j++) {
						if (AminoAcid.GetAminoAcidFromLabel(labels2[j]).Letter == aa) {
							ind2 = j;
						}
					}
					if (ind1 == -1 && ind2 == -1) {
						return -1;
					}
					if (ind1 == -1) {
						if (m != ushort.MaxValue) {
							return 2;
						}
						return -1;
					}
					if (ind2 == -1) {
						if (m != ushort.MaxValue) {
							return 1;
						}
						return -1;
					}
					SilacLabel l1 = labels1[ind1];
					SilacLabel l2 = labels2[ind2];
					if (l1 == l2) {
						if (m == ushort.MaxValue) {
							return 0;
						}
						return -1;
					}
					if (m == ushort.MaxValue) {
						return 0;
					}
					if (Tables.modifications[AminoAcid.GetMascotModificationStringForLabel(l1)].Index == m) {
						return 1;
					}
					if (Tables.modifications[AminoAcid.GetMascotModificationStringForLabel(l2)].Index == m) {
						return 2;
					}
					throw new Exception("Never get here.");
				}
				throw new Exception("Should never happen.");
			}
			return -1;
		}

		public int Charge {
			get { return charge; }
		}

		public int Id {
			get { return id; }
			set { id = value; }
		}

		public int EvidenceId {
			get { return evidenceId; }
			set { evidenceId = value; }
		}

		public double Pep {
			get { return pep; }
		}

		public double ElutionTime {
			get { return elutionTime; }
		}

		public double Score {
			get { return score; }
		}

		public double AltScore {
			get { return altScore; }
		}

		public double DeltaScoreAll {
			get { return deltaScoreAll; }
		}

		public double DeltaScoreDifferentPep {
			get { return deltaScoreDifferentPep; }
		}

		public PeptideModificationState TrueModifications {
			get { return trueModifications; }
		}

		public int ScanNumber {
			get { return scanNumber; }
		}

		public long PtmCombinatorics {
			get { return ptmScoreCounts; }
		}

		public double PtmDelta {
			get { return deltaPtmScore; }
		}

		public double PtmScore {
			get { return ptmScore; }
		}

		public double MassErrorPpm {
			get { return massErrorPpm; }
		}

		public int RawFileIndex {
			get { return rawFileIndex; }
		}

		public MsmsPeakAnnotation[] Description {
			get { return description; }
		}

		public float[] Intensities {
			get { return intensities; }
		}

		public float[] MassDiffs {
			get { return massDiffs; }
		}

		public bool Changed {
			get { return changed; }
		}

		public double Mz {
			get { return mz; }
		}

		public double MonoisotopicMz {
			get { return monoisotopicMz; }
		}

		public int SilacIndex {
			get { return silacIndex; }
		}

		public MsmsHit(BinaryReader reader) {
			rawFileIndex = reader.ReadInt32();
			scanNumber = reader.ReadInt32();
			silacIndex = reader.ReadInt32();
			byte x = reader.ReadByte();
			switch (x) {
				case 0:
					type = MascotQueryType.Silac;
					break;
				case 1:
					type = MascotQueryType.Isotope;
					break;
				case 2:
					type = MascotQueryType.Peak;
					break;
			}
			pep = reader.ReadDouble();
			score = reader.ReadDouble();
			altScore = reader.ReadDouble();
			deltaScoreAll = reader.ReadDouble();
			deltaScoreDifferentPep = reader.ReadDouble();
			labelModifications = PeptideModificationState.Read(reader);
			trueModifications = PeptideModificationState.Read(reader);
			elutionTime = reader.ReadDouble();
			id = reader.ReadInt32();
			mz = reader.ReadDouble();
			monoisotopicMz = reader.ReadDouble();
			evidenceId = reader.ReadInt32();
			ptmScore = reader.ReadDouble();
			deltaPtmScore = reader.ReadDouble();
			ptmScoreCounts = reader.ReadInt32();
			ptmScoreFinished = reader.ReadBoolean();
			massErrorPpm = reader.ReadDouble();
			changed = reader.ReadBoolean();
			int len = reader.ReadInt32();
			description = new MsmsPeakAnnotation[len];
			intensities = new float[len];
			massDiffs = new float[len];
			for (int i = 0; i < len; i++) {
				description[i] = new MsmsPeakAnnotation(reader);
				intensities[i] = reader.ReadSingle();
				massDiffs[i] = reader.ReadSingle();
			}
			charge = reader.ReadInt32();
			modProbabilities = new Dictionary<ushort, double[]>();
			len = reader.ReadInt32();
			ushort[] keys = new ushort[len];
			for (int i = 0; i < keys.Length; i++) {
				keys[i] = reader.ReadUInt16();
			}
			foreach (ushort key in keys) {
				len = reader.ReadInt32();
				double[] value = new double[len];
				for (int i = 0; i < value.Length; i++) {
					value[i] = reader.ReadDouble();
				}
				modProbabilities.Add(key, value);
			}
			modScoreDiffs = new Dictionary<ushort, double[]>();
			len = reader.ReadInt32();
			keys = new ushort[len];
			for (int i = 0; i < keys.Length; i++) {
				keys[i] = reader.ReadUInt16();
			}
			foreach (ushort key in keys) {
				len = reader.ReadInt32();
				double[] value = new double[len];
				for (int i = 0; i < value.Length; i++) {
					value[i] = reader.ReadDouble();
				}
				modScoreDiffs.Add(key, value);
			}
		}

		public void Write(BinaryWriter writer) {
			writer.Write(rawFileIndex);
			writer.Write(scanNumber);
			writer.Write(silacIndex);
			switch (type) {
				case MascotQueryType.Silac:
					writer.Write((byte) 0);
					break;
				case MascotQueryType.Isotope:
					writer.Write((byte) 1);
					break;
				case MascotQueryType.Peak:
					writer.Write((byte) 2);
					break;
			}
			writer.Write(pep);
			writer.Write(score);
			writer.Write(altScore);
			writer.Write(deltaScoreAll);
			writer.Write(deltaScoreDifferentPep);
			labelModifications.Write(writer);
			trueModifications.Write(writer);
			writer.Write(elutionTime);
			writer.Write(id);
			writer.Write(mz);
			writer.Write(monoisotopicMz);
			writer.Write(evidenceId);
			writer.Write(ptmScore);
			writer.Write(deltaPtmScore);
			writer.Write(ptmScoreCounts);
			writer.Write(ptmScoreFinished);
			writer.Write(MassErrorPpm);
			writer.Write(changed);
			writer.Write(description.Length);
			for (int i = 0; i < description.Length; i++) {
				description[i].Write(writer);
				writer.Write(intensities[i]);
				writer.Write(massDiffs[i]);
			}
			writer.Write(charge);
			ushort[] keys = DataUtil.GetKeys(modProbabilities);
			writer.Write(keys.Length);
			for (int i = 0; i < keys.Length; i++) {
				writer.Write(keys[i]);
			}
			foreach (ushort key in keys) {
				double[] value = modProbabilities[key];
				writer.Write(value.Length);
				for (int i = 0; i < value.Length; i++) {
					writer.Write(value[i]);
				}
			}
			keys = DataUtil.GetKeys(modScoreDiffs);
			writer.Write(keys.Length);
			for (int i = 0; i < keys.Length; i++) {
				writer.Write(keys[i]);
			}
			foreach (ushort key in keys) {
				double[] value = modScoreDiffs[key];
				writer.Write(value.Length);
				for (int i = 0; i < value.Length; i++) {
					writer.Write(value[i]);
				}
			}
		}

		public PeptideModificationState GetFixedAndLabelModifications(string[] fixedModifications, string sequence) {
			PeptideModificationState result = labelModifications.Clone();
			Modification[] fixedMods = Tables.ToModifications(fixedModifications);
			result.ApplyFixedModifications(fixedMods, sequence, true, true);
			return result;
		}

		public void CalcPtmScore(double ms2Tol, string ms2TolUnit, int topx, string[] fixedModifications, string sequence, double[] specMasses, float[] specIntensities) {
			List<string> fixedMods = new List<string>();
			for (int i = 0; i < fixedModifications.Length; i++) {
				fixedMods.Add(fixedModifications[i]);
			}
			ushort[] ltypes = labelModifications.ProjectToCounts().ModificationTypes;
			foreach (ushort ltype in ltypes) {
				fixedMods.Add(Tables.modificationList[ltype].Title);
			}
			PeptideModificationCounts c = trueModifications.ProjectToCounts();
			ushort[] ttypes = c.ModificationTypes;
			string[] varMods = new string[ttypes.Length];
			int[] tcounts = new int[ttypes.Length];
			for (int i = 0; i < ttypes.Length; i++) {
				varMods[i] = Tables.modificationList[ttypes[i]].Title;
				tcounts[i] = c.GetModificationCount(ttypes[i]);
			}
			PeptideModificationState modstr;
			ptmScore = Pscore.CalcPtmScore(ms2Tol, ms2TolUnit, topx, sequence, Tables.ToModifications(fixedMods.ToArray()),
			                               Tables.ToModifications(varMods), tcounts, specMasses, specIntensities,
			                               mz, charge, out ptmScoreFinished, out modstr, out ptmScoreCounts,
			                               out deltaPtmScore, out description, out intensities, out massDiffs,
			                               out modProbabilities,
			                               out modScoreDiffs, false);
			ptmScoreFinished = !ptmScoreFinished;
			changed = false;
			if (!modstr.Equals(trueModifications)) {
				changed = true;
				trueModifications = modstr;
			}
		}

		public double GetPositionalProbability(ushort mod, int pos) {
			if (!modProbabilities.ContainsKey(mod)) {
				return 0;
			}
			return modProbabilities[mod][pos];
		}

		public double GetScoreDiff(ushort mod, int pos) {
			if (!modScoreDiffs.ContainsKey(mod)) {
				return double.NaN;
			}
			return modScoreDiffs[mod][pos];
		}

		public string GetSilacStatestring(SilacType type1) {
			if (type1 != SilacType.Singlets && silacIndex != -1) {
				if (type1 == SilacType.Doublets) {
					return silacIndex == 0 ? "Light" : "Heavy";
				}
				return silacIndex == 0 ? "Light" : (silacIndex == 1 ? "Medium" : "Heavy");
			}
			return null;
		}

		public ushort[] GetAllTrueInternalModifications() {
			return trueModifications.GetAllInternalModifications();
		}

		public void Dispose() {
			labelModifications = null;
			trueModifications = null;
			description = null;
			intensities = null;
			massDiffs = null;
			if (modProbabilities != null) {
				modProbabilities.Clear();
				modProbabilities = null;
			}
			if (modScoreDiffs != null) {
				modScoreDiffs.Clear();
				modScoreDiffs = null;
			}
		}

		public double CalcMass(string sequence, string[] fixedModifications) {
			Peptide peptide = new Peptide(sequence);
			peptide.ApplyFixedModifications(Tables.ToModifications(fixedModifications));
			double result = peptide.MonoIsotopicMass;
			result += labelModifications.GetDeltaMass();
			result += trueModifications.GetDeltaMass();
			return result;
		}
	}
}