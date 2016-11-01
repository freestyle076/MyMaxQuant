/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.IO;
using MaxQuant.Mol;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public class PeptideModificationCounts {
		private int hash;
		private ushort[] modificationCounts;
		private ushort[] modificationTypes;

		public PeptideModificationCounts(ushort[] types, ushort[] counts) {
			modificationTypes = types;
			modificationCounts = counts;
		}

		public PeptideModificationCounts() {
			modificationTypes = new ushort[0];
			modificationCounts = new ushort[0];
		}

		public PeptideModificationCounts(BinaryReader reader) {
			int len = reader.ReadInt32();
			modificationTypes = new ushort[len];
			modificationCounts = new ushort[len];
			for (int i = 0; i < len; i++) {
				modificationTypes[i] = reader.ReadUInt16();
				modificationCounts[i] = reader.ReadUInt16();
			}
		}

		public double DeltaMass {
			get {
				double deltaMass = 0;
				for (int i = 0; i < modificationCounts.Length; i++) {
					Modification mod = Tables.modificationList[modificationTypes[i]];
					deltaMass += modificationCounts[i] * mod.DeltaMass;
				}
				return deltaMass;
			}
		}

		public ushort[] ModificationTypes {
			get { return modificationTypes; }
		}

		public ushort[] ModificationCounts {
			get { return modificationCounts; }
		}

		public void Write(BinaryWriter writer) {
			writer.Write(modificationTypes.Length);
			for (int i = 0; i < modificationTypes.Length; i++) {
				writer.Write(modificationTypes[i]);
				writer.Write(modificationCounts[i]);
			}
		}

		public override bool Equals(object obj) {
			if (this == obj) {
				return true;
			}
			if (obj is PeptideModificationCounts) {
				PeptideModificationCounts other = (PeptideModificationCounts) obj;
				if (other.modificationTypes.Length != modificationTypes.Length) {
					return false;
				}
				for (int i = 0; i < modificationTypes.Length; i++) {
					if (modificationTypes[i] != other.modificationTypes[i]) {
						return false;
					}
				}
				for (int i = 0; i < modificationTypes.Length; i++) {
					if (modificationCounts[i] != other.modificationCounts[i]) {
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
				for (int i = 0; i < modificationTypes.Length; i++) {
					h = 31 * h + modificationTypes[i];
				}
				for (int i = 0; i < modificationTypes.Length; i++) {
					h = 31 * h + modificationCounts[i];
				}
				hash = h;
			}
			return h;
		}

		public override string ToString() {
			if (modificationTypes.Length == 0) {
				return "Unmodified";
			}
			string[] s = new string[modificationTypes.Length];
			for (int i = 0; i < s.Length; i++) {
				if (modificationCounts[i] == 1) {
					s[i] = Tables.modificationList[modificationTypes[i]].Title;
				} else {
					s[i] = "" + modificationCounts[i] + " " + Tables.modificationList[modificationTypes[i]].Title;
				}
			}
			return StringUtil.Concat(",", s);
		}

		public ushort GetModificationCount(ushort type) {
			int index = Array.BinarySearch(modificationTypes, type);
			if (index < 0) {
				return 0;
			}
			return modificationCounts[index];
		}

		public int GetCount(ushort type) {
			int a = Array.BinarySearch(modificationTypes, type);
			if (a < 0) {
				return 0;
			}
			return modificationCounts[a];
		}

		public void Dispose() {
			modificationTypes = null;
			modificationCounts = null;
		}
	}
}