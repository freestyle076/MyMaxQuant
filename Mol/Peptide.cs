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

namespace MaxQuant.Mol {
	public class Peptide {
		private AASequence sequence;
		private double monoIsotopicMass = Double.NaN;
		private readonly List<int> proteinIndices = new List<int>(1);
		private readonly List<int> proteinOffsets = new List<int>(1);
		private readonly List<byte> residueBefore = new List<byte>(1);
		private readonly List<byte> residueAfter = new List<byte>(1);
		private PeptideModificationState fixedModifications;


		public static Peptide Read(BinaryReader reader) {
			Peptide result = new Peptide();
			result.sequence = AASequence.Read(reader);
			if (result.sequence == null) {
				return null;
			}
			result.monoIsotopicMass = reader.ReadDouble();
			int nproteins = reader.ReadInt32();
			for (int i = 0; i < nproteins; i++) {
				result.proteinIndices.Add(reader.ReadInt32());
				result.proteinOffsets.Add(reader.ReadInt32());
				result.residueBefore.Add(reader.ReadByte());
				result.residueAfter.Add(reader.ReadByte());
			}
			result.fixedModifications = PeptideModificationState.Read(reader);
			return result;
		}

		private Peptide() {
		}

		public Peptide(string sequence) {
			this.sequence = new AASequence(sequence);
			monoIsotopicMass = MonoIsotopicMass;
			fixedModifications = new PeptideModificationState(sequence.Length);
		}

		public int GetProteinIndex(int index) {
			return proteinIndices[index];
		}

		public int GetProteinOffset(int index) {
			return proteinOffsets[index];
		}

		public string Sequence {
			get { return sequence.ToString(); }
		}

		public void AddProtein(int proteinIndex, int proteinOffset, char before, char after) {
			proteinIndices.Add(proteinIndex);
			proteinOffsets.Add(proteinOffset);
			residueBefore.Add(before == ' ' ? byte.MaxValue : (byte) AminoAcid.singleLetterAas.IndexOf(before));
			residueAfter.Add(after == ' ' ? byte.MaxValue : (byte) AminoAcid.singleLetterAas.IndexOf(after));
		}

		public int GetOccurenceCount(string residues) {
			int count = 0;
			string s = Sequence;
			for (int i = 0; i < s.Length; i++) {
				if (residues.IndexOf(s[i]) != -1) {
					count++;
				}
			}
			return count;
		}

		public double MonoIsotopicMass {
			get {
				if (Double.IsNaN(monoIsotopicMass)) {
					monoIsotopicMass = MolUtil.MassNormalCTerminus + MolUtil.MassNormalNTerminus;
					string s = Sequence;
					for (int i = 0; i < s.Length; i++) {
						monoIsotopicMass += AminoAcid.aaMonoMasses[s[i]];
					}
				}
				return monoIsotopicMass;
			}
			set { monoIsotopicMass = value; }
		}

		public double ProteinCount {
			get { return proteinIndices.Count; }
		}

		public PeptideModificationState FixedModifications {
			get { return fixedModifications; }
		}

		public char GetResidueBefore(int index) {
			if (residueBefore[index] == byte.MaxValue) {
				return ' ';
			}
			return AminoAcid.singleLetterAas[residueBefore[index]];
		}

		public char GetResidueAfter(int index) {
			if (residueAfter[index] == byte.MaxValue) {
				return ' ';
			}
			return AminoAcid.singleLetterAas[residueAfter[index]];
		}

		public void Write(BinaryWriter writer) {
			sequence.Write(writer);
			writer.Write(monoIsotopicMass);
			writer.Write(proteinIndices.Count);
			for (int i = 0; i < proteinIndices.Count; i++) {
				writer.Write(proteinIndices[i]);
				writer.Write(proteinOffsets[i]);
				writer.Write(residueBefore[i]);
				writer.Write(residueAfter[i]);
			}
			fixedModifications.Write(writer);
		}

		public void ApplyFixedModifications(Modification[] modifications) {
			monoIsotopicMass += fixedModifications.ApplyFixedModifications(modifications, Sequence, IsNterm(), IsCterm());
		}

		private bool IsNterm() {
			for (int i = 0; i < ProteinCount; i++) {
				if (residueBefore[i] == ' ') {
					return true;
				}
			}
			return false;
		}

		private bool IsCterm() {
			for (int i = 0; i < ProteinCount; i++) {
				if (residueAfter[i] == ' ') {
					return true;
				}
			}
			return false;
		}

		public ModifiedPeptide SetVariableModifications(PeptideModificationState varMods, int index) {
			ModifiedPeptide result = CreateNonmodifiedVersion(index);
			result.modifications = varMods;
			if (varMods != null) {
				result.mass += varMods.GetDeltaMass();
			}
			return result;
		}

		public ModifiedPeptide[] ApplyVariableModifications(Modification[] modifications, int index) {
			ModifiedPeptide[] result = new ModifiedPeptide[] {CreateNonmodifiedVersion(index)};
			for (int i = 0; i < modifications.Length; i++) {
				result = ApplyVariableModification(result, modifications[i]);
			}
			result = FilterEqualMasses(result);
			return result;
		}

		
		public PeptideModificationState[] ApplyVariableModificationsFixedNumbers(Modification[] modifications, int[] counts,
		                                                                         out bool incomplete) {
			incomplete = false;
			PeptideModificationState[] mods = new PeptideModificationState[] {CreateNonmodifiedVersion()};
			int limit = (int)Math.Round(6000000.0 / sequence.ToString().Length / modifications.Length);
			for (int i = 0; i < modifications.Length; i++) {
				bool incompl;
				Modification mod = modifications[i];
				int n = counts[i];
				mods = ApplyVariableModificationFixedNumber(mods, mod, n, out incompl, limit);
				if (incompl) {
					incomplete = true;
				}
			}
			return mods;
		}

		private PeptideModificationState[] ApplyVariableModificationFixedNumber(IEnumerable<PeptideModificationState> mods,
		                                                                        Modification mod, int n, out bool incomplete, int limit) {
			incomplete = false;
			List<PeptideModificationState> result = new List<PeptideModificationState>();
			foreach (PeptideModificationState p in mods) {
				bool incompl;
				PeptideModificationState[] x = ApplyVariableModificationFixedNumber(p, mod, n, out incompl);
				if (incompl) {
					incomplete = true;
				}
				result.AddRange(x);
				if(result.Count > limit) {
					incomplete = true;
					break;
				}
			}
			return result.ToArray();
		}

		private static ModifiedPeptide[] FilterEqualMasses(IEnumerable<ModifiedPeptide> peptides) {
			Dictionary<double, ModifiedPeptide> x = new Dictionary<double, ModifiedPeptide>();
			foreach (ModifiedPeptide p in peptides) {
				if (!x.ContainsKey(p.mass)) {
					x.Add(p.mass, p);
				}
			}
			ModifiedPeptide[] result = new ModifiedPeptide[x.Count];
			x.Values.CopyTo(result, 0);
			return result;
		}

		private ModifiedPeptide[] ApplyVariableModification(IEnumerable<ModifiedPeptide> peptides, Modification mod) {
			List<ModifiedPeptide> result = new List<ModifiedPeptide>();
			foreach (ModifiedPeptide p in peptides) {
				ModifiedPeptide[] x = ApplyVariableModification(p, mod);
				foreach (ModifiedPeptide y in x) {
					result.Add(y);
				}
			}
			return result.ToArray();
		}

		private PeptideModificationState[] ApplyVariableModificationFixedNumber(PeptideModificationState peptide,
		                                                                        Modification mod, int n,
		                                                                        out bool incomplete) {
			incomplete = false;
			if (n == 0) {
				return new PeptideModificationState[0];
			}
			if (mod.GetPosition() == ModificationPosition.anyNterm) {
				if (n > 1) {
					return new PeptideModificationState[0];
				}
				if (peptide.GetNTermModification() == ushort.MaxValue) {
					PeptideModificationState q = peptide.Clone();
					q.SetNTermModification(mod.Index);
					return new PeptideModificationState[] {q};
				}
				return new PeptideModificationState[0];
			}
			if (mod.GetPosition() == ModificationPosition.anyCterm) {
				if (n > 1) {
					return new PeptideModificationState[0];
				}
				if (peptide.GetCTermModification() == ushort.MaxValue) {
					PeptideModificationState q = peptide.Clone();
					q.SetCTermModification(mod.Index);
					return new PeptideModificationState[] {q};
				}
				return new PeptideModificationState[0];
			}
			if (mod.GetPosition() == ModificationPosition.proteinNterm) {
				if (n > 1) {
					return new PeptideModificationState[0];
				}
				if (peptide.GetNTermModification() == ushort.MaxValue) {
					PeptideModificationState q = peptide.Clone();
					q.SetNTermModification(mod.Index);
					return new PeptideModificationState[] {q};
				}
				return new PeptideModificationState[0];
			}
			if (mod.GetPosition() == ModificationPosition.proteinCterm) {
				if (n > 1) {
					return new PeptideModificationState[0];
				}
				if (peptide.GetCTermModification() == ushort.MaxValue) {
					PeptideModificationState q = peptide.Clone();
					q.SetCTermModification(mod.Index);
					return new PeptideModificationState[] {q};
				}
				return new PeptideModificationState[0];
			}
			string s = Sequence;
			switch (mod.GetTermType(0)) {
				case ModificationSiteType.aa: {
					HashSet<int> ind = new HashSet<int>();
					for (int j = 0; j < s.Length; j++) {
						for (int i = 0; i < mod.AaCount; i++) {
							if (s[j] == mod.GetAaAt(i) && peptide.GetModificationAt(j) == ushort.MaxValue) {
								ind.Add(j);
							}
						}
					}
					int[] indices = ind.ToArray();
					if (indices.Length < n) {
						return new PeptideModificationState[0];
					}
					List<PeptideModificationState> result = new List<PeptideModificationState>();
					bool incompl;
					int[][] comb = NumUtil.GetCombinations(indices.Length, n, 10000, out incompl);
					if (incompl) {
						incomplete = incompl;
					}
					for (int j = 0; j < comb.Length; j++) {
						int[] iind = ArrayUtil.SubArray(indices, comb[j]);
						PeptideModificationState q = peptide.Clone();
						for (int k = 0; k < iind.Length; k++) {
							q.SetModificationAt(iind[k], mod.Index);
						}
						result.Add(q);
					}
					return result.ToArray();
				}
				case ModificationSiteType.nterm: {
					if (n > 1) {
						return new PeptideModificationState[0];
					}
					for (int i = 0; i < mod.AaCount; i++) {
						List<PeptideModificationState> toBeAdded = new List<PeptideModificationState>();
						if (s[0] == mod.GetAaAt(i) && peptide.GetNTermModification() == ushort.MaxValue) {
							PeptideModificationState q = peptide.Clone();
							q.SetNTermModification(mod.Index);
							toBeAdded.Add(q);
						}
						if (toBeAdded.Count > 0) {
							return toBeAdded.ToArray();
						}
					}
					return new PeptideModificationState[0];
				}
				case ModificationSiteType.cterm: {
					if (n > 1) {
						return new PeptideModificationState[0];
					}
					for (int i = 0; i < mod.AaCount; i++) {
						List<PeptideModificationState> toBeAdded = new List<PeptideModificationState>();
						if (s[s.Length - 1] == mod.GetAaAt(i) && peptide.GetCTermModification() == ushort.MaxValue) {
							PeptideModificationState q = peptide.Clone();
							q.SetCTermModification(mod.Index);
							toBeAdded.Add(q);
						}
						if (toBeAdded.Count > 0) {
							return toBeAdded.ToArray();
						}
					}
					return new PeptideModificationState[0];
				}
			}
			throw new Exception("Never get here.");
		}

		private ModifiedPeptide[] ApplyVariableModification(ModifiedPeptide peptide, Modification mod) {
			List<ModifiedPeptide> result = new List<ModifiedPeptide>();
			result.Add(peptide);
			string s = Sequence;
			for (int i = 0; i < mod.AaCount; i++) {
				switch (mod.GetTermType(i)) {
					case ModificationSiteType.aa: {
						List<ModifiedPeptide> toBeAdded = new List<ModifiedPeptide>();
						foreach (ModifiedPeptide w in result) {
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
							result.Add(a);
						}
						break;
					}
					case ModificationSiteType.nterm: {
						List<ModifiedPeptide> toBeAdded = new List<ModifiedPeptide>();
						foreach (ModifiedPeptide w in result) {
							if (s[0] == mod.GetAaAt(i) && w.modifications.GetNTermModification() == ushort.MaxValue) {
								ModifiedPeptide q = w.Clone();
								q.mass += mod.DeltaMass;
								q.modifications.SetNTermModification(mod.Index);
								toBeAdded.Add(q);
							}
						}
						foreach (ModifiedPeptide a in toBeAdded) {
							result.Add(a);
						}
						break;
					}
					case ModificationSiteType.cterm: {
						List<ModifiedPeptide> toBeAdded = new List<ModifiedPeptide>();
						foreach (ModifiedPeptide w in result) {
							if (s[s.Length - 1] == mod.GetAaAt(i) && w.modifications.GetCTermModification() == ushort.MaxValue) {
								ModifiedPeptide q = w.Clone();
								q.mass += mod.DeltaMass;
								q.modifications.SetCTermModification(mod.Index);
								toBeAdded.Add(q);
							}
						}
						foreach (ModifiedPeptide a in toBeAdded) {
							result.Add(a);
						}
						break;
					}
				}
			}
			if (mod.GetPosition() == ModificationPosition.anyNterm) {
				List<ModifiedPeptide> toBeAdded = new List<ModifiedPeptide>();
				foreach (ModifiedPeptide w in result) {
					if (w.modifications.GetNTermModification() == ushort.MaxValue) {
						ModifiedPeptide q = w.Clone();
						q.mass += mod.DeltaMass;
						q.modifications.SetNTermModification(mod.Index);
						toBeAdded.Add(q);
					}
				}
				foreach (ModifiedPeptide a in toBeAdded) {
					result.Add(a);
				}
			}
			if (mod.GetPosition() == ModificationPosition.anyCterm) {
				List<ModifiedPeptide> toBeAdded = new List<ModifiedPeptide>();
				foreach (ModifiedPeptide w in result) {
					if (w.modifications.GetCTermModification() == ushort.MaxValue) {
						ModifiedPeptide q = w.Clone();
						q.mass += mod.DeltaMass;
						q.modifications.SetCTermModification(mod.Index);
						toBeAdded.Add(q);
					}
				}
				foreach (ModifiedPeptide a in toBeAdded) {
					result.Add(a);
				}
			}
			if (mod.GetPosition() == ModificationPosition.proteinNterm && IsNterm()) {
				List<ModifiedPeptide> toBeAdded = new List<ModifiedPeptide>();
				foreach (ModifiedPeptide w in result) {
					if (w.modifications.GetNTermModification() == ushort.MaxValue) {
						ModifiedPeptide q = w.Clone();
						q.mass += mod.DeltaMass;
						q.modifications.SetNTermModification(mod.Index);
						toBeAdded.Add(q);
					}
				}
				foreach (ModifiedPeptide a in toBeAdded) {
					result.Add(a);
				}
			}
			if (mod.GetPosition() == ModificationPosition.proteinCterm && IsCterm()) {
				List<ModifiedPeptide> toBeAdded = new List<ModifiedPeptide>();
				foreach (ModifiedPeptide w in result) {
					if (w.modifications.GetCTermModification() == ushort.MaxValue) {
						ModifiedPeptide q = w.Clone();
						q.mass += mod.DeltaMass;
						q.modifications.SetCTermModification(mod.Index);
						toBeAdded.Add(q);
					}
				}
				foreach (ModifiedPeptide a in toBeAdded) {
					result.Add(a);
				}
			}
			return result.ToArray();
		}

		private PeptideModificationState CreateNonmodifiedVersion() {
			return new PeptideModificationState(Sequence.Length);
		}

		private ModifiedPeptide CreateNonmodifiedVersion(int index) {
			return new ModifiedPeptide(this, index);
		}
	}
}