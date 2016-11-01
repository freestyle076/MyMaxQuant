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
	[Serializable]
	public class MascotPeptide {
		private readonly double mass;//peptide mass, not m/z
		private float deltaMass;
		private string sequence;
		private PeptideModificationState modifications;
		private float mascotScore;
		private float altScore = float.NaN;
		private double pep = 1;
		private int[] proteinIndex;

		public MascotPeptide(string[] peptideInfo, int[] proteinIndex, string[] modstrings) {
			mass = Double.Parse(peptideInfo[1]);
			sequence = peptideInfo[4].ToUpper();
			modifications = new PeptideModificationState(peptideInfo[6], modstrings);
			mascotScore = (float) Double.Parse(peptideInfo[7]);
			this.proteinIndex = proteinIndex;
			deltaMass = (float)Double.Parse(peptideInfo[2]);
		}

		public void Correct(string prefix, IProteinSet proteinSet) {
			int count = 0;
			for (int i = 0; i < proteinIndex.Length; i++) {
				if (T07SearchEngineEnhancement.IsReverseProtein(proteinSet.GetName(proteinIndex[i]), prefix)) {
					count++;
				}
			}
			if (count == 0 || count == proteinIndex.Length) {
				return;
			}
			List<int> result = new List<int>();
			for (int i = 0; i < proteinIndex.Length; i++) {
				if (!proteinSet.GetName(proteinIndex[i]).StartsWith(prefix)) {
					result.Add(proteinIndex[i]);
				}
			}
			proteinIndex = result.ToArray();
		}

		public PeptideModificationState GetTrueModifications(HashSet<string> labelModificationSet) {
			return Modifications.GetTrueModifications(labelModificationSet);
		}

		public string GetModificationClassification(HashSet<string> labelModificationSet) {
			PeptideModificationState trueMods = GetTrueModifications(labelModificationSet);
			ushort[] x = trueMods.GetModifications();
			int index = -1;
			int count = 0;
			for(int i = 0; i < x.Length; i++) {
				if(x[i] != ushort.MaxValue) {
					count++;
					index = i;
				}
			}
			if(count ==0) {
				return "unmodified";
			}
			if(count > 1) {
				return "multiple";
			}
			return "" + x[index] + ";" + sequence[index];
		}

		public MascotPeptide(BinaryReader reader) {
			mass = reader.ReadDouble();
			deltaMass = reader.ReadSingle();
			sequence = reader.ReadString();
			modifications = PeptideModificationState.Read(reader);
			mascotScore = reader.ReadSingle();
			altScore = reader.ReadSingle();
			pep = reader.ReadDouble();
			int len = reader.ReadInt32();
			proteinIndex = new int[len];
			for (int i = 0; i < len; i++) {
				proteinIndex[i] = reader.ReadInt32();
			}
		}

		public string Sequence {
			get { return sequence; }
		}

		public int[] ProteinIndex {
			get { return proteinIndex; }
		}

		public double MascotScore {
			get { return mascotScore; }
			set { mascotScore = (float) value; }
		}

		public bool Unique {
			get { return proteinIndex.Length == 1; }
		}

		public PeptideModificationState Modifications {
			get { return modifications; }
			set { modifications = value; }
		}

		public void Write(BinaryWriter writer) {
			writer.Write(mass);
			writer.Write(deltaMass);
			writer.Write(sequence);
			if (modifications == null) {
				modifications = new PeptideModificationState(sequence.Length);
			}
			modifications.Write(writer);
			writer.Write(mascotScore);
			writer.Write(altScore);
			writer.Write(pep);
			int len = proteinIndex.Length;
			writer.Write(len);
			for (int i = 0; i < len; i++) {
				writer.Write(proteinIndex[i]);
			}
		}

		public bool FitsAaCounts(char[] aa, int[] count) {
			for (int i = 0; i < aa.Length; i++) {
				if (!FitsAaCount(aa[i], count[i])) {
					return false;
				}
			}
			return true;
		}

		public bool FitsAaCount(char aa, int count) {
			int c = 0;
			string s = Sequence;
			for (int i = 0; i < s.Length; i++) {
				if (s[i] == aa) {
					c++;
				}
			}
			return count == c;
		}

		public int GetAaCount(char aa) {
			int c = 0;
			string s = Sequence;
			for (int i = 0; i < s.Length; i++) {
				if (s[i] == aa) {
					c++;
				}
			}
			return c;
		}

		public bool HasOnlyReverseHits(string reverseStr, IProteinSet proteinSet) {
			for (int i = 0; i < proteinIndex.Length; i++) {
				if (!T07SearchEngineEnhancement.IsReverseProtein(proteinSet.GetName(proteinIndex[i]), reverseStr)) {
					return false;
				}
			}
			return true;
		}

		public string[] GetProteinIds(IProteinSet proteinSet) {
			string[] result = new string[proteinIndex.Length];
			for (int i = 0; i < result.Length; i++) {
				result[i] = proteinSet.GetName(proteinIndex[i]);
			}
			return result;
		}

		public double MassErrorPpm {
			get { return deltaMass / mass * 1e6; }
		}

		public double Mass {
			get { return mass; }
		}

		public double DeltaMass {
			get { return deltaMass; }
			set { deltaMass = (float) value; }
		}

		public double Pep {
			get { return pep; }
			set { pep = value; }
		}

		public double AltScore {
			get { return altScore; }
			set { altScore = (float) value; }
		}

		public int GetSilacIndex(char[] aas, SilacType type, SilacLabel[] labels1, SilacLabel[] labels2) {
			foreach (char aa in aas) {
				int ind = GetSilacIndex(aa, type, labels1, labels2);
				if (ind != -1) {
					return ind;
				}
			}
			return -1;
		}

		private int GetSilacIndex(char aa, SilacType type, SilacLabel[] labels1, SilacLabel[] labels2) {
			if (type == SilacType.Singlets) {
				return 0;
			}
			int index = -1;
			string seq = sequence;
			for (int i = 0; i < seq.Length; i++) {
				if (seq[i] == aa) {
					index = i;
					break;
				}
			}
			if (index == -1) {
				return -1;
			}
			Modification[] mods = T07SearchEngineEnhancement.GetModifiedVersions(aa, type, labels1, labels2);
			ushort mod = modifications.GetModificationAt(index);
			if (mod == ushort.MaxValue) {
				for (int i = 1; i < mods.Length; i++) {
					if (mods[i] == null) {
						return -1;
					}
				}
				return 0;
			}
			if (type == SilacType.Doublets) {
				if (mods[1] == null) {
					return -1;
				}
				if (mods[1].Index == mod) {
					return 1;
				}
				return -1;
			}
			if (mods[1] == mods[2]) {
				return -1;
			}
			if (mods[1] != null && mods[1].Index == mod) {
				return 1;
			}
			if (mods[2] != null && mods[2].Index == mod) {
				return 2;
			}
			return -1;
		}

		public bool HasConsistentLabelModifications(char[] aas, HashSet<string> labelModificationSet) {
			foreach (char aa in aas) {
				if (!HasConsistentLabelModifications(aa, labelModificationSet)) {
					return false;
				}
			}
			return true;
		}

		private bool HasConsistentLabelModifications(char aa, HashSet<string> labelModificationSet) {
			List<ushort> mods = new List<ushort>();
			string seq = sequence;
			for (int i = 0; i < seq.Length; i++) {
				if (seq[i] == aa) {
					mods.Add(modifications.GetModificationAt(i));
				}
			}
			if (mods.Count < 2) {
				return true;
			}
			int labelIndex = -1;
			for (int i = 0; i < mods.Count; i++) {
				if (mods[i] != ushort.MaxValue && Tables.modificationList[mods[i]].IsLabelModification(labelModificationSet)) {
					labelIndex = i;
					break;
				}
			}
			if (labelIndex == -1) {
				return true;
			}
			ushort labelMod = mods[labelIndex];
			for (int i = 0; i < mods.Count; i++) {
				if (mods[i] != labelMod) {
					return false;
				}
			}
			return true;
		}

		public void Dispose() {
			if(sequence != null) {
				//sequence.Dispose();
				sequence = null;						
			}
			if (modifications != null) {
				modifications.Dispose();
				modifications = null;
			}
			proteinIndex = null;
		}

		public bool ValidModifications() {
			PeptideModificationCounts counts = modifications.ProjectToCounts();
			ushort[] c = counts.ModificationCounts;
			foreach (ushort cc in c) {
				if (cc > 7) {
					return false;
				}
			}
			return true;
		}

		public static double CalcMonoisotopicMass(string sequence, PeptideModificationState varMod, Modification[] fixMod) {
			double result = AminoAcid.CalcMonoisotopicMass(sequence) + varMod.GetDeltaMass();
			foreach (char c in sequence) {
				result += GetDeltaFixed(c, fixMod);
			}
			return result;
		}

		private static double GetDeltaFixed(char c, IEnumerable<Modification> mod) {
			foreach (Modification m in mod) {
				for(int i = 0; i < m.AaCount; i++) {
					if(c == m.GetAaAt(i)) {
						return m.DeltaMass;
					}
				}
			}
			return 0;
		}

		public void RecalcMassDeviation(ushort[] mods, IPeakList peakList, int scanNumber, SilacType type, double isotopeCorrelationThreshold) {
			Modification[] fixedMods = ArrayUtil.SubArray(Tables.modificationList, mods);
			double m = CalcMonoisotopicMass(sequence, modifications, fixedMods);
			int charge;
			double mz;
			double mass1;
			T06MsmsPreparation.GetPrecursorInfo(out charge, out mz, out mass1, peakList, scanNumber, type, isotopeCorrelationThreshold, double.NaN);
			if(!double.IsNaN(mass1) && mass1 > 0){
				float x = (float) (mass1 - m);
				if(!float.IsNaN(x) && Math.Abs(x/m)<1e-5 ){
					deltaMass = x;
				}
			}	
		}
	}
}