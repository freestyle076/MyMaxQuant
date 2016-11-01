/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.IO;
using System.Text;

namespace MaxQuant.Mol {
	[Serializable]
	public class AASequence : IComparable {
		private const int aasPerLong = 12;
		private const int basis = 32;

		private ulong[] sequence;
		private int len;
		private int hash;

		public static AASequence Read(BinaryReader reader) {
			try {
				AASequence result = new AASequence();
				result.len = reader.ReadInt32();
				result.sequence = new ulong[reader.ReadInt32()];
				for (int i = 0; i < result.sequence.Length; i++) {
					result.sequence[i] = reader.ReadUInt64();
				}
				return result;
			}catch(EndOfStreamException) {
				return null;
			}
		}

		public AASequence(){}

		public AASequence(string seq) {
			len = seq.Length;
			int a = len / aasPerLong;
			int b = len % aasPerLong;
			int n = (b == 0) ? a : a + 1;
			sequence = new ulong[n];
			for (int i = 0; i < a; i++) {
				sequence[i] = Encode(seq.Substring(aasPerLong * i, aasPerLong));
			}
			if (b > 0) {
				sequence[a] = Encode(seq.Substring(aasPerLong * a, len - aasPerLong * a));
			}
		}

		private static ulong Encode(string seq) {
			ulong result = 0;
			for (int i = 0; i < seq.Length; i++) {
				char c = seq[i];
				int x = AminoAcid.singleLetterAas.IndexOf(c);
				if (x == -1) {
					x = 20;
				}
				result = (result << 5) + (ulong) x;
			}
			return result;
		}

		private static string Decode(ulong a, int n) {
			char[] result = new char[n];
			for (int i = 0; i < n; i++) {
				int x = (int) (a % basis);
				a = a >> 5;
				result[n - 1 - i] = (x == 20) ? 'X' : AminoAcid.singleLetterAas[x];
			}
			return new string(result);
		}

		public override string ToString() {
			StringBuilder b = new StringBuilder();
			for (int i = 0; i < sequence.Length - 1; i++) {
				b.Append(Decode(sequence[i], aasPerLong));
			}
			int l = len - (sequence.Length - 1) * aasPerLong;
			b.Append(Decode(sequence[sequence.Length - 1], l));
			return b.ToString();
		}

		public override int GetHashCode() {
			int h = hash;
			if (h == 0) {
				for (int i = 0; i < sequence.Length; i++) {
					h = 31 * h + (int) sequence[i];
				}
				hash = h;
			}
			return h;
		}

		public override bool Equals(object obj) {
			if (this == obj) {
				return true;
			}
			if (obj is AASequence) {
				AASequence other = (AASequence) obj;
				if (other.len != len) {
					return false;
				}
				for (int i = 0; i < sequence.Length; i++) {
					if (sequence[i] != other.sequence[i]) {
						return false;
					}
				}
				return true;
			}
			return false;
		}

		public void Write(BinaryWriter writer) {
			writer.Write(len);
			writer.Write(sequence.Length);
			foreach (ulong u in sequence) {
				writer.Write(u);
			}
		}

		public int CompareTo(object obj) {
			AASequence other = (AASequence) obj;
			int min = Math.Min(sequence.Length, other.sequence.Length);
			for (int i = 0; i < min; i++) {
				if (sequence[i] != other.sequence[i]) {
					if (sequence[i] < other.sequence[i]) {
						return -1;
					}
					return 1;
				}
			}
			if (sequence.Length != other.sequence.Length) {
				if (sequence.Length < other.sequence.Length) {
					return -1;
				}
				return 1;
			}
			if (len != other.len) {
				if (len < other.len) {
					return -1;
				}
				return 1;
			}
			return 0;
		}

		public void Dispose() {
			sequence = null;
		}
	}
}