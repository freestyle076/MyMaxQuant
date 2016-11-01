/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Mol {
	public class Modification {
		private double deltaMass = double.NaN;
		private ushort index;
		private double[] neutralLoss;
		private ModificationPosition position = ModificationPosition.anywhere;
		private char[] site;
		private ModificationSiteType[] termType;
		private string title;

		public int AaCount {
			get { return site.Length; }
		}

		public double DeltaMass {
			get { return deltaMass; }
			set { deltaMass = value; }
		}

		public string Abbreviation {
			get { return title.Substring(0, 2).ToLower(); }
		}

		public bool HasNeutralLoss {
			get { return neutralLoss.Length > 0; }
		}

		public bool IsPhosphorylation {
			get { return title.ToLower().Contains("phospho"); }
		}

		public bool IsInternal {
			get {
				if (termType.Length == 0) {
					return false;
				}
				return termType[0] == ModificationSiteType.aa;
			}
		}

		public ushort Index {
			get { return index; }
			set { index = value; }
		}

		public string Title {
			get { return title; }
			set { title = value; }
		}

		public char[] Site {
			set { site = value; }
		}

		public ModificationPosition GetPosition() {
			return position;
		}

		public void SetPosition(ModificationPosition pos) {
			position = pos;
		}

		public ModificationSiteType GetTermType(int i) {
			return termType[i];
		}

		public override bool Equals(object obj) {
			if (this == obj) {
				return true;
			}
			if (obj is Modification) {
				Modification other = (Modification) obj;
				return (other.index != index);
			}
			return false;
		}

		public override int GetHashCode() {
			return index;
		}

		public bool IsLabelModification(HashSet<string> labelModificationSet) {
			return labelModificationSet.Contains(title);
		}

		public bool HasAA(char aa) {
			foreach (char x in site) {
				if (x == aa) {
					return true;
				}
			}
			return false;
		}

		public char GetAaAt(int j) {
			return site[j];
		}

		public void SetTermType(ModificationSiteType[] termTypes) {
			termType = termTypes;
		}

		public void SetNeutralLoss(double[] x) {
			neutralLoss = x;
		}

		public double[] GetNeutralLoss() {
			return neutralLoss;
		}
	}
}