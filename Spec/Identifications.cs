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
	public class Identifications {
		private string filePath;
		private ushort[][] fixedModifications;
		private MascotPeptide[][] peptides;
		private int[] scanNumbers;
		private MascotQueryType type;

		public Identifications(string filePath) {
			this.filePath = filePath;
		}

		public Identifications(MascotQuery[] queries, string filePath, MascotQueryType type) {
			this.filePath = filePath;
			this.type = type;
			scanNumbers = new int[queries.Length];
			peptides = new MascotPeptide[queries.Length][];
			fixedModifications = new ushort[queries.Length][];
			for (int i = 0; i < scanNumbers.Length; i++) {
				scanNumbers[i] = queries[i].ScanNum;
			}
			int[] o = ArrayUtil.Order(scanNumbers);
			scanNumbers = ArrayUtil.SubArray(scanNumbers, o);
			queries = ArrayUtil.SubArray(queries, o);
			for (int i = 0; i < queries.Length; i++) {
				peptides[i] = queries[i].GetPeptides();
				fixedModifications[i] = queries[i].FixedMods;
				double[] scores = new double[peptides[i].Length];
				for (int j = 0; j < scores.Length; j++) {
					scores[j] = -peptides[i][j].MascotScore;
				}
				int[] o1 = ArrayUtil.Order(scores);
				peptides[i] = ArrayUtil.SubArray(peptides[i], o1);
			}
			List<int> valids = new List<int>();
			for (int i = 0; i < scanNumbers.Length; i++) {
				if (peptides[i].Length > 0) {
					valids.Add(i);
				}
			}
			int[] oo = valids.ToArray();
			scanNumbers = ArrayUtil.SubArray(scanNumbers, oo);
			peptides = ArrayUtil.SubArray(peptides, oo);
			fixedModifications = ArrayUtil.SubArray(fixedModifications, oo);
		}

		public int Count {
			get {
				if (scanNumbers == null) {
					Read();
				}
				return scanNumbers.Length;
			}
		}

		public int[] ScanNumbers {
			get {
				if (scanNumbers == null) {
					Read();
				}
				return scanNumbers;
			}
		}

		public MascotPeptide[][] Peptides {
			get {
				if (peptides == null) {
					Read();
				}
				return peptides;
			}
		}

		public ushort[][] FixedModifications {
			get { return fixedModifications; }
		}

		public MascotPeptide[] GetPeptidesAt(int index) {
			if (peptides == null) {
				Read();
			}
			return peptides[index];
		}

		public double GetHighestMascotScore(int index) {
			if (peptides == null) {
				Read();
			}
			return peptides[index][0].MascotScore;
		}

		public double GetHighestAltScore(int index) {
			if (peptides == null) {
				Read();
			}
			return peptides[index][0].AltScore;
		}

		public string GetBestSequence(int index) {
			if (peptides == null) {
				Read();
			}
			return peptides[index][0].Sequence;
		}

		public bool IsHighestScoringCorrect(int index, string revstring, IProteinSet proteinSet) {
			MascotPeptide p = GetPeptidesAt(index)[0];
			return !p.HasOnlyReverseHits(revstring, proteinSet);
		}

		public string GetModificationClassification(int index, HashSet<string> labelModSet) {
			MascotPeptide p = GetPeptidesAt(index)[0];
			return p.GetModificationClassification(labelModSet);
		}

		public MascotPeptide GetTopPeptideByScanNumber(int scanNumber) {
			MascotPeptide[] peps = GetPeptidesByScanNumber(scanNumber);
			if (peps == null) {
				return null;
			}
			return peps[0];
		}

		public MascotPeptide[] GetPeptidesByScanNumber(int scanNumber) {
			int ind = GetIndexByScanNumber(scanNumber);
			if (ind < 0) {
				return null;
			}
			return peptides[ind];
		}

		public int GetIndexByScanNumber(int scanNumber) {
			if (peptides == null) {
				Read();
			}
			int ind = Array.BinarySearch(scanNumbers, scanNumber);
			if (ind < 0) {
				return -1;
			}
			return ind;
		}

		public void Read() {
			if (!File.Exists(filePath)) {
				return;
			}
			using (BinaryReader reader = FileUtil.GetBinaryReader(filePath)) {
				int len = reader.ReadInt32();
				scanNumbers = new int[len];
				peptides = new MascotPeptide[len][];
				fixedModifications = new ushort[len][];
				for (int i = 0; i < scanNumbers.Length; i++) {
					scanNumbers[i] = reader.ReadInt32();
					int l = reader.ReadInt32();
					peptides[i] = new MascotPeptide[l];
					for (int j = 0; j < peptides[i].Length; j++) {
						peptides[i][j] = new MascotPeptide(reader);
					}
					l = reader.ReadInt32();
					fixedModifications[i] = new ushort[l];
					for (int j = 0; j < fixedModifications[i].Length; j++) {
						fixedModifications[i][j] = reader.ReadUInt16();
					}
				}
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
				reader.Close();
			}
		}

		public void Write() {
			if (File.Exists(filePath)) {
				File.Delete(filePath);
			}
			using (BinaryWriter writer = FileUtil.GetBinaryWriter(filePath)) {
				writer.Write(scanNumbers.Length);
				for (int i = 0; i < scanNumbers.Length; i++) {
					writer.Write(scanNumbers[i]);
					writer.Write(peptides[i].Length);
					for (int j = 0; j < peptides[i].Length; j++) {
						peptides[i][j].Write(writer);
					}
					writer.Write(fixedModifications[i].Length);
					for (int j = 0; j < fixedModifications[i].Length; j++) {
						writer.Write(fixedModifications[i][j]);
					}
				}
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
				writer.Flush();
				writer.Close();
			}
		}

		public void SetData(int[] scanNumbersIn, MascotPeptide[][] peptidesIn, ushort[][] fixedModificationsIn) {
			scanNumbers = scanNumbersIn;
			peptides = peptidesIn;
			fixedModifications = fixedModificationsIn;
		}

		public void Extract(int[] indices) {
			scanNumbers = ArrayUtil.SubArray(scanNumbers, indices);
			peptides = ArrayUtil.SubArray(peptides, indices);
			fixedModifications = ArrayUtil.SubArray(fixedModifications, indices);
		}

		public void FillProteinToPepTable(Dictionary<string, HashSet<string>> protIdToPepSeqs, IProteinSet proteinSet) {
			if (peptides == null) {
				Read();
			}
			for (int i = 0; i < peptides.Length; i++) {
				MascotPeptide p = peptides[i][0];
				int[] proteinIndex = p.ProteinIndex;
				if (!T07SearchEngineEnhancement.ValidPeptide(p.Sequence)) {
					continue;
				}
				foreach (int pi in proteinIndex) {
					string protId = proteinSet.GetName(pi);
					if (!protIdToPepSeqs.ContainsKey(protId)) {
						protIdToPepSeqs.Add(protId, new HashSet<string>());
					}
					string key = p.Sequence;
					if (!protIdToPepSeqs[protId].Contains(key)) {
						protIdToPepSeqs[protId].Add(key);
					}
				}
			}
		}

		public bool Exists() {
			return File.Exists(filePath);
		}

		public void ProcessPeptides(int fileIndex, Dictionary<string, int> proteinIdToGroupIndex, IPeakList peakList,
		                            MsmsData msmsData, IIdentifiedPeptide[] identifiedPeptides, string[] peptideSequences,
		                            ReQuantitationResult reQuantitationResult, bool reQuantify,
		                            HashSet<string> labelModificationSet, IProteinSet proteinSet,
		                            SilacType silacType, SilacLabel[] labels1, SilacLabel[] labels2, double ms2Tol,
		                            string ms2TolUnit, int topx, string[] fixedMods) {
			if (peptides == null) {
				Read();
			}
			double[] monoIsoMz = peakList.MS2MonoisotopicMz;
			for (int i = 0; i < peptides.Length; i++) {
				MascotPeptide[] p = peptides[i];
				int scanNumber = scanNumbers[i];
				int ms2ind = peakList.GetMs2IndexFromScanNumber(scanNumber);
				double mz = peakList.GetMs2Mz(ms2ind);
				double monotopicMz = monoIsoMz[ms2ind];
				double time = peakList.GetMs2Rt(ms2ind);
				int silacId = -1;
				int isotopeId = -1;
				int silacIndex = -1;
				SilacCluster silacCluster = null;
				IsotopeCluster isotopeCluster = null;
				if (type == MascotQueryType.Silac) {
					int[] silacInfo = peakList.GetSilacInfoForMsmsScanNumber(scanNumber);
					silacId = silacInfo[0];
					silacIndex = silacInfo[1];
					silacCluster = peakList.GetSilacCluster(silacId);
				} else if (type == MascotQueryType.Isotope) {
					isotopeId = peakList.GetIsotopeIndexForMsmsScanNumber(scanNumber);
					isotopeCluster = peakList.GetIsotopeCluster(isotopeId);
				}
				int index = Array.BinarySearch(peptideSequences, p[0].Sequence);
				if (index < 0) {
					continue;
				}
				HashSet<int> tmpGroupInds = new HashSet<int>();
				foreach (int pi in p[0].ProteinIndex) {
					string protId = proteinSet.GetName(pi);
					if (!proteinIdToGroupIndex.ContainsKey(protId)) {
						continue;
					}
					int groupInd = proteinIdToGroupIndex[protId];
					if (!tmpGroupInds.Contains(groupInd)) {
						tmpGroupInds.Add(groupInd);
					}
				}
				double[] specMasses;
				float[] specIntensities;
				bool uniqueProtein = (p[0].ProteinIndex.Length == 1);
				bool uniqueGroup = (tmpGroupInds.Count == 1);
				msmsData.GetSpectrumFromScanNumber(scanNumber, out specMasses, out specIntensities);
				identifiedPeptides[index].AddMascotPeptideHit(p, scanNumber, fileIndex, type, silacId, silacIndex, isotopeId,
				                                              silacCluster, isotopeCluster, time, peakList, mz, monotopicMz,
				                                              fixedModifications[i], specMasses, specIntensities,
				                                              reQuantitationResult, reQuantify, labelModificationSet, silacType,
				                                              labels1, labels2, ms2Tol, ms2TolUnit, topx, fixedMods);
				identifiedPeptides[index].UniqueProtein = uniqueProtein;
				identifiedPeptides[index].UniqueGroup = uniqueGroup;
			}
		}

		public void Dispose() {
			filePath = null;
			scanNumbers = null;
			if (peptides != null) {
				for (int i = 0; i < peptides.Length; i++) {
					for (int j = 0; j < peptides[i].Length; j++) {
						peptides[i][j] = null;
					}
					peptides[i] = null;
				}
			}
			peptides = null;
			if (fixedModifications != null) {
				for (int i = 0; i < fixedModifications.Length; i++) {
					fixedModifications[i] = null;
				}
			}
			fixedModifications = null;
		}

		public void DeepDispose() {
			filePath = null;
			scanNumbers = null;
			if (peptides != null) {
				for (int i = 0; i < peptides.Length; i++) {
					if (peptides[i] != null) {
						for (int j = 0; j < peptides[i].Length; j++) {
							if (peptides[i][j] != null) {
								peptides[i][j].Dispose();
								peptides[i][j] = null;
							}
						}
						peptides[i] = null;
					}
				}
			}
			peptides = null;
			if (fixedModifications != null) {
				for (int i = 0; i < fixedModifications.Length; i++) {
					fixedModifications[i] = null;
				}
				fixedModifications = null;
			}
		}

		public Modification[] GetFixedModifications(int ind) {
			Modification[] result = new Modification[fixedModifications[ind].Length];
			for (int i = 0; i < fixedModifications[ind].Length; i++) {
				result[i] = Tables.modificationList[fixedModifications[ind][i]];
			}
			return result;
		}

		public void SetPeptidesAt(int i, MascotPeptide[] mascotPeptides) {
			peptides[i] = mascotPeptides;
		}

		public int GetScanNumberAt(int i) {
			return scanNumbers[i];
		}

		public void Merge(Identifications other) {
			if (type != other.type) {
				throw new Exception("Different types.");
			}
			scanNumbers = ArrayUtil.Concat(scanNumbers, other.scanNumbers);
			peptides = ArrayUtil.Concat(peptides, other.peptides);
			fixedModifications = ArrayUtil.Concat(fixedModifications, other.fixedModifications);
			int[] o = ArrayUtil.Order(scanNumbers);
			scanNumbers = ArrayUtil.SubArray(scanNumbers, o);
			peptides = ArrayUtil.SubArray(peptides, o);
			fixedModifications = ArrayUtil.SubArray(fixedModifications, o);
		}

		public void AddSequences(HashSet<string> x) {
			if (peptides == null) {
				Read();
			}
			for (int i = 0; i < peptides.Length; i++) {
				MascotPeptide[] p = peptides[i];
				x.Add(p[0].Sequence);
			}
		}

		public bool Contains(int scanNum) {
			return Array.BinarySearch(ScanNumbers, scanNum) >= 0;
		}
	}
}