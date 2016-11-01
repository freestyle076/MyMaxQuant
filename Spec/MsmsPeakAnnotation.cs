/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.IO;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public class MsmsPeakAnnotation {
		private readonly MsmsSeriesType type;
		private readonly int index;
		private readonly int charge;
		private readonly bool neutralLoss;

		/// <summary>
		/// Cache the hash code for the <c>PeptideModificationState</c>
		/// </summary>
		private int hash;

		public MsmsPeakAnnotation(MsmsSeriesType type, int index, int charge, bool neutralLoss) {
			this.type = type;
			this.index = index;
			this.charge = charge;
			this.neutralLoss = neutralLoss;
		}

		public MsmsPeakAnnotation(BinaryReader reader) {
			int a = reader.ReadInt32();
			switch(a) {
				case 0:
					type = MsmsSeriesType.A;
					break;
				case 1:
					type = MsmsSeriesType.B;
					break;
				case 2:
					type = MsmsSeriesType.C;
					break;
				case 3:
					type = MsmsSeriesType.X;
					break;
				case 4:
					type = MsmsSeriesType.Y;
					break;
				case 5:
					type = MsmsSeriesType.Z;
					break;
			}
			index = reader.ReadInt32();
			charge = reader.ReadInt32();
			neutralLoss = reader.ReadBoolean();
		}

		public MsmsSeriesType Type {
			get { return type; }
		}

		public bool NeutralLoss {
			get { return neutralLoss; }
		}

		public int Index {
			get { return index; }
		}

		public int Charge {
			get { return charge; }
		}

		public bool IsNTerminal() {
			return type == MsmsSeriesType.A || type == MsmsSeriesType.B || type == MsmsSeriesType.C;
		}

		public bool IsCTerminal() {
			return type == MsmsSeriesType.X || type == MsmsSeriesType.Y || type == MsmsSeriesType.Z;
		}

		public void Write(BinaryWriter writer) {
			switch(type) {
				case MsmsSeriesType.A:
					writer.Write(0);
					break;
				case MsmsSeriesType.B:
					writer.Write(1);
					break;
				case MsmsSeriesType.C:
					writer.Write(2);
					break;
				case MsmsSeriesType.X:
					writer.Write(3);
					break;
				case MsmsSeriesType.Y:
					writer.Write(4);
					break;
				case MsmsSeriesType.Z:
					writer.Write(5);
					break;
			}
			writer.Write(index);
			writer.Write(charge);
			writer.Write(neutralLoss);
		}

		public static string[] ToStrings2(MsmsPeakAnnotation[] annot) {
			string[] result = new string[annot.Length];
			for (int i = 0; i < result.Length; i++) {
				result[i] = annot[i].ToString2();
			}
			return result;
		}

		public static string[] ToStrings(MsmsPeakAnnotation[] annot) {
			string[] result = new string[annot.Length];
			for (int i = 0; i < result.Length; i++) {
				result[i] = annot[i].ToString();
			}
			return result;
		}

		public override string ToString() {
			string result = null;
			switch (type) {
				case MsmsSeriesType.A:
					result = "a";
					break;
				case MsmsSeriesType.B:
					result = "b";
					break;
				case MsmsSeriesType.C:
					result = "c";
					break;
				case MsmsSeriesType.X:
					result = "x";
					break;
				case MsmsSeriesType.Y:
					result = "y";
					break;
				case MsmsSeriesType.Z:
					result = "z";
					break;
			}
			result += index;
			if (charge > 1) {
				result += "(" + charge + "+)";
			}
			if (neutralLoss) {
				result += "*";
			}
			return result;
		}

		public string ToString2() {
			string result = null;
			switch (type) {
				case MsmsSeriesType.A:
					result = "a";
					break;
				case MsmsSeriesType.B:
					result = "b";
					break;
				case MsmsSeriesType.C:
					result = "c";
					break;
				case MsmsSeriesType.X:
					result = "x";
					break;
				case MsmsSeriesType.Y:
					result = "y";
					break;
				case MsmsSeriesType.Z:
					result = "z";
					break;
			}
			result += StringUtil.ToSubscript(index, false);
			if (charge > 1) {
				result += StringUtil.ToSuperscript(charge, false) + '\u207A';
			}
			if (neutralLoss) {
				result += "*";
			}
			return result;
		}

		public override bool Equals(object anObject) {
			if (this == anObject) {
				return true;
			}
			if (anObject is MsmsPeakAnnotation) {
				MsmsPeakAnnotation other = (MsmsPeakAnnotation)anObject;
				if (other.type != type) {
					return false;
				}
				if (other.index != index) {
					return false;
				}
				if (other.charge != charge) {
					return false;
				}
				if (other.neutralLoss != neutralLoss) {
					return false;
				}
				return true;
			}
			return false;
		}

		public override int GetHashCode() {
			int h = hash;
			if (h == 0) {
				h = 31 * h + index;
				h = 31 * h + charge;
				h = 31 * h + (int)type;
				h = 31 * h + (neutralLoss ? 1 : 0);
				hash = h;
			}
			return h;
		}
	}
}