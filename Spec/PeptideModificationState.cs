/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.IO;
using MaxQuant.Mol;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	[Serializable]
	public class PeptideModificationState {

		/// <summary>
		/// Cache the hash code for the <c>PeptideModificationState</c>
		/// </summary>
		private int hash;

		private ushort[] modifications;
		private ushort nTermModification = ushort.MaxValue;
		private ushort cTermModification = ushort.MaxValue;

		private PeptideModificationState() {
		}

		public PeptideModificationState(int len) {
			modifications = new ushort[len];
			for (int i = 0; i < len; i++) {
				modifications[i] = ushort.MaxValue;
			}
		}

		public PeptideModificationState(string mascotString, string[] modStrings)
			: this(mascotString.Length - 2) {
			if (mascotString[0] != '0') {
				int ind = Int32.Parse("" + mascotString[0]) - 1;
				nTermModification = Tables.modifications[modStrings[ind]].Index;
			}
			if (mascotString[mascotString.Length - 1] != '0') {
				int ind = Int32.Parse("" + mascotString[mascotString.Length - 1]) - 1;
				cTermModification = Tables.modifications[modStrings[ind]].Index;
			}
			for (int i = 0; i < modifications.Length; i++) {
				if (mascotString[i + 1] != '0') {
					int ind = Int32.Parse("" + mascotString[i + 1]) - 1;
					modifications[i] = Tables.modifications[modStrings[ind]].Index;
				}
			}
		}

		public int Length {
			get { return modifications.Length; }
		}

		public int Count {
			get {
				int c = 0;
				if (nTermModification != ushort.MaxValue) {
					c++;
				}
				if (cTermModification != ushort.MaxValue) {
					c++;
				}
				foreach (ushort m in modifications) {
					if (m != ushort.MaxValue) {
						c++;
					}
				}
				return c;
			}
		}

		public ushort GetNTermModification() {
			return nTermModification;
		}

		public void SetNTermModification(ushort value) {
			nTermModification = value;
		}

		public ushort GetCTermModification() {
			return cTermModification;
		}

		public void SetCTermModification(ushort value) {
			cTermModification = value;
		}

		public ushort[] GetModifications() {
			return modifications;
		}

		public ushort GetModificationAt(int index) {
			return modifications[index];
		}

		public void SetModificationAt(int index, ushort value) {
			modifications[index] = value;
		}

		public static PeptideModificationState Read(BinaryReader reader) {
			PeptideModificationState result = new PeptideModificationState();
			result.cTermModification = reader.ReadUInt16();
			result.nTermModification = reader.ReadUInt16();
			int len = reader.ReadInt32();
			result.modifications = new ushort[len];
			for (int i = 0; i < result.modifications.Length; i++) {
				result.modifications[i] = reader.ReadUInt16();
			}
			return result;
		}

		public void Write(BinaryWriter writer) {
			writer.Write(cTermModification);
			writer.Write(nTermModification);
			writer.Write(modifications.Length);
			for (int i = 0; i < modifications.Length; i++) {
				writer.Write(modifications[i]);
			}
		}

		public long SizeEstimate() {
			return 8 + 2 * modifications.Length;
		}

		public PeptideModificationState Clone() {
			PeptideModificationState p = new PeptideModificationState();
			p.modifications = (ushort[]) modifications.Clone();
			p.nTermModification = nTermModification;
			p.cTermModification = cTermModification;
			return p;
		}

		public PeptideModificationState GetFreshCopy(int len) {
			PeptideModificationState result = new PeptideModificationState();
			result.SetNTermModification(ushort.MaxValue);
			result.SetCTermModification(ushort.MaxValue);
			result.modifications = new ushort[len];
			for (int i = 0; i < len; i++) {
				result.modifications[i] = ushort.MaxValue;
			}
			return result;
		}

		public PeptideModificationCounts ProjectToCounts() {
			return ProjectToCounts(modifications, nTermModification, cTermModification);
		}

		public static PeptideModificationCounts ProjectToCounts(ushort[] modifications, ushort nTermModification,
		                                                        ushort cTermModification) {
			ushort[] c = new ushort[Tables.modificationList.Length];
			for (int i = 0; i < modifications.Length; i++) {
				if (modifications[i] != ushort.MaxValue) {
					c[modifications[i]]++;
				}
			}
			if (nTermModification != ushort.MaxValue) {
				c[nTermModification]++;
			}
			if (cTermModification != ushort.MaxValue) {
				c[cTermModification]++;
			}
			List<ushort> types = new List<ushort>();
			List<ushort> counts = new List<ushort>();
			for (int i = 0; i < c.Length; i++) {
				if (c[i] > 0) {
					types.Add((ushort) i);
					counts.Add(c[i]);
				}
			}
			return new PeptideModificationCounts(types.ToArray(), counts.ToArray());
		}

		public ushort[] GetAllInternalModifications() {
			List<ushort> result = new List<ushort>();
			ushort[] b = GetModifications();
			for (int i = 0; i < b.Length; i++) {
				if (b[i] != ushort.MaxValue) {
					result.Add(b[i]);
				}
			}
			return result.ToArray();
		}

		public override bool Equals(object anObject) {
			if (this == anObject) {
				return true;
			}
			if (anObject is PeptideModificationState) {
				PeptideModificationState other = (PeptideModificationState) anObject;
				if (other.GetNTermModification() != GetNTermModification()) {
					return false;
				}
				if (other.GetCTermModification() != GetCTermModification()) {
					return false;
				}
				if (other.Length != Length) {
					return false;
				}
				for (int i = 0; i < Length; i++) {
					if (other.GetModificationAt(i) != GetModificationAt(i)) {
						return false;
					}
				}
				return true;
			}
			return false;
		}

		public override int GetHashCode() {
			int h = hash;
			if (h == 0) {
				for (int i = 0; i < Length; i++) {
					h = 31 * h + GetModificationAt(i);
				}
				h = 31 * h + GetNTermModification();
				h = 31 * h + GetCTermModification();
				hash = h;
			}
			return h;
		}

		public double ApplyFixedModifications(Modification[] mods, string sequence, bool isNterm, bool isCterm) {
			double deltaMass = 0;
			foreach (Modification mod in mods) {
				deltaMass += ApplyFixedModification(mod, sequence, isNterm, isCterm);
			}
			return deltaMass;
		}

		private double ApplyFixedModification(Modification mod, string sequence, bool isNterm, bool isCterm) {
			double deltaMass = 0;
			for (int i = 0; i < mod.AaCount; i++) {
				switch (mod.GetTermType(i)) {
					case ModificationSiteType.aa:
						for (int j = 0; j < sequence.Length; j++) {
							if (sequence[j] == mod.GetAaAt(i)) {
								if (GetModificationAt(j) != ushort.MaxValue) {
									throw new Exception("Conflicting fixed modifications.");
								}
								SetModificationAt(j, mod.Index);
								deltaMass += mod.DeltaMass;
							}
						}
						break;
					case ModificationSiteType.nterm:
						if (sequence[0] == mod.GetAaAt(i)) {
							if (GetNTermModification() != ushort.MaxValue) {
								throw new Exception("Conflicting fixed modifications.");
							}
							SetNTermModification(mod.Index);
							deltaMass += mod.DeltaMass;
						}
						break;
					case ModificationSiteType.cterm:
						if (sequence[sequence.Length - 1] == mod.GetAaAt(i)) {
							if (GetCTermModification() != ushort.MaxValue) {
								throw new Exception("Conflicting fixed modifications.");
							}
							SetCTermModification(mod.Index);
							deltaMass += mod.DeltaMass;
						}
						break;
				}
			}
			if (mod.GetPosition() == ModificationPosition.anyNterm) {
				if (GetNTermModification() != ushort.MaxValue) {
					throw new Exception("Conflicting fixed modifications.");
				}
				SetNTermModification(mod.Index);
				deltaMass += mod.DeltaMass;
			}
			if (mod.GetPosition() == ModificationPosition.anyCterm) {
				if (GetCTermModification() != ushort.MaxValue) {
					throw new Exception("Conflicting fixed modifications.");
				}
				SetCTermModification(mod.Index);
				deltaMass += mod.DeltaMass;
			}
			if (mod.GetPosition() == ModificationPosition.proteinNterm && isNterm) {
				if (GetNTermModification() != ushort.MaxValue) {
					throw new Exception("Conflicting fixed modifications.");
				}
				SetNTermModification(mod.Index);
				deltaMass += mod.DeltaMass;
			}
			if (mod.GetPosition() == ModificationPosition.proteinCterm && isCterm) {
				if (GetCTermModification() != ushort.MaxValue) {
					throw new Exception("Conflicting fixed modifications.");
				}
				SetCTermModification(mod.Index);
				deltaMass += mod.DeltaMass;
			}
			return deltaMass;
		}

		public PeptideModificationState GetTrueModifications(HashSet<string> labelModificationSet) {
			PeptideModificationState result = GetFreshCopy(Length);
			result.SetNTermModification(GetNTermModification());
			result.SetCTermModification(GetCTermModification());
			for (int i = 0; i < Length; i++) {
				ushort m = GetModificationAt(i);
				if (m == ushort.MaxValue || labelModificationSet.Contains(Tables.modificationList[m].Title)) {
					result.SetModificationAt(i, ushort.MaxValue);
				} else {
					result.SetModificationAt(i, GetModificationAt(i));
				}
			}
			return result;
		}

		public PeptideModificationState GetLabelModifications(ushort[] fixedMods, string sequence,
		                                                      HashSet<string> labelModificationSet) {
			PeptideModificationState result = GetFreshCopy(Length);
			for (int i = 0; i < Length; i++) {
				ushort m = GetModificationAt(i);
				if (m != ushort.MaxValue && labelModificationSet.Contains(Tables.modificationList[m].Title)) {
					result.SetModificationAt(i, GetModificationAt(i));
				} else {
					result.SetModificationAt(i, ushort.MaxValue);
				}
			}
			foreach (ushort fixMod in fixedMods) {
				if (labelModificationSet.Contains(Tables.modificationList[fixMod].Title)) {
					Modification mod = Tables.modificationList[fixMod];
					for (int j = 0; j < mod.AaCount; j++) {
						char c = mod.GetAaAt(j);
						for (int i = 0; i < Length; i++) {
							if (sequence[i] == c) {
								result.SetModificationAt(i, mod.Index);
							}
						}
					}
				}
			}
			return result;
		}

		public double GetDeltaMass() {
			double result = 0;
			if (nTermModification != ushort.MaxValue) {
				result += Tables.modificationList[nTermModification].DeltaMass;
			}
			if (cTermModification != ushort.MaxValue) {
				result += Tables.modificationList[cTermModification].DeltaMass;
			}
			foreach (ushort m in modifications) {
				if (m != ushort.MaxValue) {
					result += Tables.modificationList[m].DeltaMass;
				}
			}
			return result;
		}

		public PeptideModificationState RemoveLabelModifications(HashSet<string> labelModificationSet) {
			PeptideModificationState result = Clone();
			if (nTermModification != ushort.MaxValue &&
			    Tables.modificationList[nTermModification].IsLabelModification(labelModificationSet)) {
				result.nTermModification = ushort.MaxValue;
			}
			if (cTermModification != ushort.MaxValue &&
			    Tables.modificationList[cTermModification].IsLabelModification(labelModificationSet)) {
				result.cTermModification = ushort.MaxValue;
			}
			for (int i = 0; i < modifications.Length; i++) {
				if (modifications[i] != ushort.MaxValue &&
				    Tables.modificationList[modifications[i]].IsLabelModification(labelModificationSet)) {
					result.modifications[i] = ushort.MaxValue;
				}
			}
			return result;
		}

		public PeptideModificationState Revert() {
			PeptideModificationState result = Clone();
			for (int i = 0; i < modifications.Length; i++) {
				result.modifications[i] = modifications[modifications.Length - 1 - i];
			}
			return result;
		}

		public void Dispose() {
			modifications = null;
		}
	}
}