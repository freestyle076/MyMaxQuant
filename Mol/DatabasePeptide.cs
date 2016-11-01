/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.IO;
using MaxQuant.Mol;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Mol {
	public class DatabasePeptide {
		private byte length;
		private double monoisotopicMass = Double.NaN;
		private readonly List<int> proteinIndices = new List<int>(1);
		private readonly List<ushort> proteinOffsets = new List<ushort>(1);
		private PeptideModificationCounts fixedModifications = new PeptideModificationCounts();

		public static DatabasePeptide Read(BinaryReader reader) {
			try {
				DatabasePeptide result = new DatabasePeptide();
				result.length = reader.ReadByte();
				result.monoisotopicMass = reader.ReadDouble();
				int nproteins = reader.ReadInt32();
				for (int i = 0; i < nproteins; i++) {
					result.proteinIndices.Add(reader.ReadInt32());
					result.proteinOffsets.Add(reader.ReadUInt16());
				}
				result.fixedModifications = new PeptideModificationCounts(reader);
				return result;
			} catch (EndOfStreamException) {
				return null;
			}
		}

		private DatabasePeptide() {
		}

		public DatabasePeptide(string sequence, int proteinIndex, int proteinOffset) {
			length = (byte) sequence.Length;
			monoisotopicMass = GetMonoIsotopicMass(sequence);
			AddProtein(proteinIndex, proteinOffset);
		}

		public int GetProteinIndex(int index) {
			return proteinIndices[index];
		}

		public int[] GetProteinIndices() {
			return proteinIndices.ToArray();
		}

		public int GetProteinOffset(int index) {
			return proteinOffsets[index];
		}

		public string GetSequence(IProteinSet proteinSet) {
			string seq = proteinSet.Get(proteinIndices[0]).Sequence;
			return seq.Substring(proteinOffsets[0], length);
		}

		public void AddProtein(int proteinIndex, int proteinOffset) {
			proteinIndices.Add(proteinIndex);
			proteinOffsets.Add((ushort) proteinOffset);
		}

		public int GetOccurenceCount(string residues, IProteinSet proteinSet) {
			int count = 0;
			string s = GetSequence(proteinSet);
			for (int i = 0; i < s.Length; i++) {
				if (residues.IndexOf(s[i]) != -1) {
					count++;
				}
			}
			return count;
		}

		public double GetMonoIsotopicMass(string sequence) {
			if (Double.IsNaN(monoisotopicMass)) {
				monoisotopicMass = MolUtil.MassNormalCTerminus + MolUtil.MassNormalNTerminus;
				for (int i = 0; i < sequence.Length; i++) {
					monoisotopicMass += AminoAcid.aaMonoMasses[sequence[i]];
				}
			}
			return monoisotopicMass;
		}

		public void SetMonoisotopicMass(double monoMass) {
			monoisotopicMass = monoMass;
		}

		public double ProteinCount {
			get { return proteinIndices.Count; }
		}

		public int Length {
			get { return length; }
		}

		public PeptideModificationCounts GetModifications(double modifiedMass, ushort[] varMods, double[] deltaM, ushort[] max) {
			double dm = modifiedMass - monoisotopicMass;
			ushort[] counts = new ushort[max.Length];
			bool finished = false;
			while (!finished) {
				double x = 0;
				for (int i = 0; i < counts.Length; i++) {
					x += counts[i] * deltaM[i];
				}
				if (Math.Abs(x - dm) < 1e-6) {
					return new PeptideModificationCounts(varMods, counts);
				}
				finished = IncrementCounter(counts, max);
			}
			throw new Exception("Unable to reconstruct modifications.");
		}

		private static bool IncrementCounter(ushort[] counts, ushort[] max) {
			int index = 0;
			while (counts[index] == max[index]) {
				index++;
				if (index == counts.Length) {
					return true;
				}
			}
			counts[index]++;
			for (int i = 0; i < index; i++) {
				counts[i] = 0;
			}
			return false;
		}

		public char GetResidueBefore(int index, IProteinSet proteinSet) {
			int pos = proteinOffsets[index] - 1;
			if (pos < 0) {
				return ' ';
			}
			return proteinSet.Get(proteinIndices[index]).Sequence[pos];
		}

		public char GetResidueAfter(int index, IProteinSet proteinSet) {
			int pos = proteinOffsets[index] + length;
			string seq = proteinSet.Get(proteinIndices[index]).Sequence;
			if (pos >= seq.Length) {
				return ' ';
			}
			return seq[pos];
		}

		public void Write(BinaryWriter writer) {
			writer.Write(length);
			writer.Write(monoisotopicMass);
			writer.Write(proteinIndices.Count);
			for (int i = 0; i < proteinIndices.Count; i++) {
				writer.Write(proteinIndices[i]);
				writer.Write(proteinOffsets[i]);
			}
			fixedModifications.Write(writer);
		}

		public void ApplyFixedModifications(Modification[] modifications, IProteinSet proteinSet, string sequence) {
			PeptideModificationCounts modCounts;
			monoisotopicMass +=
				ApplyFixedModifications(modifications, sequence, IsNterm(), IsCterm(proteinSet), length, out modCounts);
			fixedModifications = modCounts;
		}

		public static double ApplyFixedModifications(Modification[] mods, string sequence, bool isNterm, bool isCterm, int len,
		                                             out PeptideModificationCounts modCounts) {
			ushort[] modifications = new ushort[len];
			for (int i = 0; i < len; i++) {
				modifications[i] = ushort.MaxValue;
			}
			ushort nTermModification = ushort.MaxValue;
			ushort cTermModification = ushort.MaxValue;
			double deltaMass = 0;
			foreach (Modification mod in mods) {
				deltaMass += ApplyFixedModification(mod, sequence, isNterm, isCterm, modifications,
				                                    ref nTermModification, ref cTermModification);
			}
			modCounts = PeptideModificationState.ProjectToCounts(modifications, nTermModification, cTermModification);
			return deltaMass;
		}

		private static double ApplyFixedModification(Modification mod, string sequence, bool isNterm, bool isCterm,
		                                             ushort[] modifications,
		                                             ref ushort nTermModification, ref ushort cTermModification) {
			double deltaMass = 0;
			for (int i = 0; i < mod.AaCount; i++) {
				switch (mod.GetTermType(i)) {
					case ModificationSiteType.aa:
						for (int j = 0; j < sequence.Length; j++) {
							if (sequence[j] == mod.GetAaAt(i)) {
								if (modifications[j] != ushort.MaxValue) {
									throw new Exception("Conflicting fixed modifications.");
								}
								modifications[j] = mod.Index;
								deltaMass += mod.DeltaMass;
							}
						}
						break;
					case ModificationSiteType.nterm:
						if (sequence[0] == mod.GetAaAt(i)) {
							if (nTermModification != ushort.MaxValue) {
								throw new Exception("Conflicting fixed modifications.");
							}
							nTermModification = mod.Index;
							deltaMass += mod.DeltaMass;
						}
						break;
					case ModificationSiteType.cterm:
						if (sequence[sequence.Length - 1] == mod.GetAaAt(i)) {
							if (cTermModification != ushort.MaxValue) {
								throw new Exception("Conflicting fixed modifications.");
							}
							cTermModification = mod.Index;
							deltaMass += mod.DeltaMass;
						}
						break;
				}
			}
			if (mod.GetPosition() == ModificationPosition.anyNterm) {
				if (nTermModification != ushort.MaxValue) {
					throw new Exception("Conflicting fixed modifications.");
				}
				nTermModification = mod.Index;
				deltaMass += mod.DeltaMass;
			}
			if (mod.GetPosition() == ModificationPosition.anyCterm) {
				if (cTermModification != ushort.MaxValue) {
					throw new Exception("Conflicting fixed modifications.");
				}
				cTermModification = mod.Index;
				deltaMass += mod.DeltaMass;
			}
			if (mod.GetPosition() == ModificationPosition.proteinNterm && isNterm) {
				if (nTermModification != ushort.MaxValue) {
					throw new Exception("Conflicting fixed modifications.");
				}
				nTermModification = mod.Index;
				deltaMass += mod.DeltaMass;
			}
			if (mod.GetPosition() == ModificationPosition.proteinCterm && isCterm) {
				if (cTermModification != ushort.MaxValue) {
					throw new Exception("Conflicting fixed modifications.");
				}
				cTermModification = mod.Index;
				deltaMass += mod.DeltaMass;
			}
			return deltaMass;
		}
		
		private bool IsNterm() {
			for (int i = 0; i < ProteinCount; i++) {
				if (proteinOffsets[i] == 0) {
					return true;
				}
			}
			return false;
		}

		private bool IsCterm(IProteinSet proteinSet) {
			for (int i = 0; i < ProteinCount; i++) {
				int pos = proteinOffsets[i] + length;
				int len = proteinSet.GetLength(proteinIndices[i]);
				if (pos >= len) {
					return true;
				}
			}
			return false;
		}

		public DatabaseModifiedPeptide[] ApplyVariableModifications(Modification[] modifications, Modification[][] lMods,
		                                                            int index, IProteinSet proteinSet) {
			string sequence = GetSequence(proteinSet);
			ModifiedPeptide[] result = new ModifiedPeptide[] { CreateNonmodifiedVersion(index, sequence) };
			result = ApplyLabelModifications(result, lMods, sequence);
			for (int i = 0; i < modifications.Length; i++) {
				result = ApplyVariableModification(result, modifications[i], sequence, proteinSet);
			}
			result = FilterEqualMods(result);
			return ConvertToDatabasePeptides(result);
		}

		private static DatabaseModifiedPeptide[] ConvertToDatabasePeptides(ModifiedPeptide[] w) {
			List<DatabaseModifiedPeptide> result = new List<DatabaseModifiedPeptide>();
			for (int i = 0; i < w.Length; i++) {
				result.Add(new DatabaseModifiedPeptide(w[i].peptideIndex, w[i].mass));
			}
			return result.ToArray();
		}

		private static ModifiedPeptide[] FilterEqualMods(IEnumerable<ModifiedPeptide> peptides) {
			Dictionary<PeptideModificationCounts, ModifiedPeptide> x =
				new Dictionary<PeptideModificationCounts, ModifiedPeptide>();
			foreach (ModifiedPeptide p in peptides) {
				PeptideModificationCounts mod = p.modifications.ProjectToCounts();
				if (!x.ContainsKey(mod)) {
					x.Add(mod, p);
				}
			}
			return DataUtil.GetValues(x);
		}

		private ModifiedPeptide[] ApplyVariableModification(IEnumerable<ModifiedPeptide> peptides, Modification mod,
		                                                    string sequence, IProteinSet proteinSet) {
			List<ModifiedPeptide> result = new List<ModifiedPeptide>();
			foreach (ModifiedPeptide p in peptides) {
				ModifiedPeptide[] x = ApplyVariableModification(p, mod, sequence, proteinSet);
				foreach (ModifiedPeptide y in x) {
					result.Add(y);
				}
			}
			//if (result.Count > 500000) {
			//    return FilterEqualMods(result.ToArray());
			//} else {
			return result.ToArray();
			//}
		}

		private readonly List<ModifiedPeptide> help = new List<ModifiedPeptide>();
		private readonly List<ModifiedPeptide> toBeAdded = new List<ModifiedPeptide>();

		private ModifiedPeptide[] ApplyVariableModification(ModifiedPeptide peptide, Modification mod, string s, IProteinSet proteinSet) {
			help.Clear();
			help.Add(peptide);
			for (int i = 0; i < mod.AaCount; i++) {
				switch (mod.GetTermType(i)) {
					case ModificationSiteType.aa: {
						toBeAdded.Clear();
						foreach (ModifiedPeptide w in help) {
							List<int> indices = new List<int>();
							for (int j = 0; j < s.Length; j++) {
								if (s[j] == mod.GetAaAt(i) && w.modifications.GetModificationAt(j) == ushort.MaxValue) {
									indices.Add(j);
								}
							}
							for (int j = 1; j <= indices.Count; j++) {
								ModifiedPeptide q = w.Clone();
								q.mass += j * mod.DeltaMass;
								for (int k = 0; k < j; k++) {
									q.modifications.SetModificationAt(indices[k], mod.Index);
								}
								toBeAdded.Add(q);
							}
						}
						foreach (ModifiedPeptide a in toBeAdded) {
							help.Add(a);
						}
						break;
					}
					case ModificationSiteType.nterm: {
						toBeAdded.Clear();
						foreach (ModifiedPeptide w in help) {
							if (s[0] == mod.GetAaAt(i) && w.modifications.GetNTermModification() == ushort.MaxValue) {
								ModifiedPeptide q = w.Clone();
								q.mass += mod.DeltaMass;
								q.modifications.SetNTermModification(mod.Index);
								toBeAdded.Add(q);
							}
						}
						foreach (ModifiedPeptide a in toBeAdded) {
							help.Add(a);
						}
						break;
					}
					case ModificationSiteType.cterm: {
						toBeAdded.Clear();
						foreach (ModifiedPeptide w in help) {
							if (s[s.Length - 1] == mod.GetAaAt(i) && w.modifications.GetCTermModification() == ushort.MaxValue) {
								ModifiedPeptide q = w.Clone();
								q.mass += mod.DeltaMass;
								q.modifications.SetCTermModification(mod.Index);
								toBeAdded.Add(q);
							}
						}
						foreach (ModifiedPeptide a in toBeAdded) {
							help.Add(a);
						}
						break;
					}
				}
			}
			if (mod.GetPosition() == ModificationPosition.anyNterm) {
				toBeAdded.Clear();
				foreach (ModifiedPeptide w in help) {
					if (w.modifications.GetNTermModification() == ushort.MaxValue) {
						ModifiedPeptide q = w.Clone();
						q.mass += mod.DeltaMass;
						q.modifications.SetNTermModification(mod.Index);
						toBeAdded.Add(q);
					}
				}
				foreach (ModifiedPeptide a in toBeAdded) {
					help.Add(a);
				}
			}
			if (mod.GetPosition() == ModificationPosition.anyCterm) {
				toBeAdded.Clear();
				foreach (ModifiedPeptide w in help) {
					if (w.modifications.GetCTermModification() == ushort.MaxValue) {
						ModifiedPeptide q = w.Clone();
						q.mass += mod.DeltaMass;
						q.modifications.SetCTermModification(mod.Index);
						toBeAdded.Add(q);
					}
				}
				foreach (ModifiedPeptide a in toBeAdded) {
					help.Add(a);
				}
			}
			if (mod.GetPosition() == ModificationPosition.proteinNterm && IsNterm()) {
				toBeAdded.Clear();
				foreach (ModifiedPeptide w in help) {
					if (w.modifications.GetNTermModification() == ushort.MaxValue) {
						ModifiedPeptide q = w.Clone();
						q.mass += mod.DeltaMass;
						q.modifications.SetNTermModification(mod.Index);
						toBeAdded.Add(q);
					}
				}
				foreach (ModifiedPeptide a in toBeAdded) {
					help.Add(a);
				}
			}
			if (mod.GetPosition() == ModificationPosition.proteinCterm && IsCterm(proteinSet)) {
				toBeAdded.Clear();
				foreach (ModifiedPeptide w in help) {
					if (w.modifications.GetCTermModification() == ushort.MaxValue) {
						ModifiedPeptide q = w.Clone();
						q.mass += mod.DeltaMass;
						q.modifications.SetCTermModification(mod.Index);
						toBeAdded.Add(q);
					}
				}
				foreach (ModifiedPeptide a in toBeAdded) {
					help.Add(a);
				}
			}
			return help.ToArray();
		}

		private static ModifiedPeptide[] ApplyLabelModifications(IEnumerable<ModifiedPeptide> peptides, Modification[][] mod,
		                                                         string sequence) {
			List<ModifiedPeptide> result = new List<ModifiedPeptide>();
			foreach (ModifiedPeptide p in peptides) {
				ModifiedPeptide[] x = ApplyLabelModifications(p, mod, sequence);
				result.AddRange(x);
			}
			return result.ToArray();
		}

		private static ModifiedPeptide[] ApplyLabelModifications(ModifiedPeptide peptide, Modification[][] mod,
		                                                         string sequence) {
			ModifiedPeptide[] result = new ModifiedPeptide[mod.Length + 1];
			result[0] = peptide;
			for (int i = 0; i < mod.Length; i++) {
				result[i + 1] = ApplyLabelModifications(peptide, mod[i], sequence);
			}
			return result;
		}

		private static ModifiedPeptide ApplyLabelModifications(ModifiedPeptide peptide, Modification[] mod,
		                                                       string s) {
			ModifiedPeptide result = peptide.Clone();
			for (int i = 0; i < mod.Length; i++) {
				for (int j = 0; j < mod[i].AaCount; j++) {
					for (int k = 0; k < s.Length; k++) {
						if (s[k] == mod[i].GetAaAt(j) && result.modifications.GetModificationAt(k) == ushort.MaxValue) {
							result.mass += mod[i].DeltaMass;
							result.modifications.SetModificationAt(k, mod[i].Index);
						}
					}
				}
			}
			return result;
		}

		private ModifiedPeptide CreateNonmodifiedVersion(int index, string sequence) {
			return new ModifiedPeptide(this, index, sequence);
		}
	}
}