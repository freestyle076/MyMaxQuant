/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Collections.Generic;
using System.IO;
using System.Text;
using MaxQuant.Mol;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Tasks {
	public static class T06MsmsPreparation {
		public static void ParseMsm(string msmName, string[] silacMsmNames, string isoMsmName, string peakMsmName,
		                            string polyMsmName, IPeakList peakList, out int silacCount, out int isoCount,
		                            out int peakCount,
		                            out int polyCount, SilacType silacType, double isotopeCorrelationThreshold) {
			HashSet<int> scanNumbers = new HashSet<int>();
			StreamReader reader = new StreamReader(msmName);
			StreamWriter[] silacWriters = new StreamWriter[silacMsmNames.GetLength(0)];
			for (int i = 0; i < silacWriters.GetLength(0); i++) {
				silacWriters[i] = new StreamWriter(silacMsmNames[i]);
			}
			StreamWriter isoWriter = new StreamWriter(isoMsmName);
			StreamWriter peakWriter = new StreamWriter(peakMsmName);
			StreamWriter polyWriter = new StreamWriter(polyMsmName);
			string line;
			List<string> buffer = new List<string>();
			string title = null;
			string mass = null;
			silacCount = 0;
			isoCount = 0;
			peakCount = 0;
			polyCount = 0;
			while ((line = reader.ReadLine()) != null) {
				if (line.IndexOf("BEGIN IONS") != -1) {
					buffer = new List<string>();
				} else if (line.IndexOf("PEPMASS") != -1) {
					mass = line;
				} else if (line.IndexOf("CHARGE") != -1) {
				} else if (line.IndexOf("TITLE") != -1) {
					title = line;
				} else if (line.IndexOf("END IONS") != -1) {
					ProcessMsm(buffer.ToArray(), title, mass, scanNumbers, silacWriters, isoWriter, peakWriter, polyWriter,
					           peakList, isotopeCorrelationThreshold, silacType,
					           ref silacCount, ref isoCount, ref peakCount, ref polyCount);
				} else {
					if (line.Length > 0) {
						buffer.Add(line);
					}
				}
			}
			reader.Close();
			for (int i = 0; i < silacWriters.GetLength(0); i++) {
				silacWriters[i].Close();
			}
			isoWriter.Close();
			peakWriter.Close();
			polyWriter.Close();
		}

		private static void ProcessMsm(string[] lines, string title, string origMass,
		                               HashSet<int> scanNumbers, StreamWriter[] silacWriters,
		                               TextWriter isoWriter, TextWriter peakWriter, TextWriter polyWriter,
		                               IPeakList peakList, double isotopeCorrelationThreshold, SilacType type,
		                               ref int silacCount, ref int isoCount, ref int peakCount, ref int polyCount) {
			int scanNumber = ExtractScanNumber(title);
			if (scanNumbers.Contains(scanNumber)) {
				return;
			}
			int[] silacInfo = peakList.GetSilacInfoForMsmsScanNumber(scanNumber);
			if (silacInfo != null) {
				scanNumbers.Add(scanNumber);
				int silacId = silacInfo[0];
				int silacIndex = silacInfo[1];
				int ch = peakList.GetSilacCharge(silacId);
				double mass = peakList.GetSilacMass(silacId, silacIndex, type);
				double mz = mass / ch + MolUtil.MassProton;
				silacWriters[silacIndex].WriteLine("BEGIN IONS");
				silacWriters[silacIndex].WriteLine("PEPMASS=" + mz);
				silacWriters[silacIndex].WriteLine("CHARGE=" + ch + "+");
				silacWriters[silacIndex].WriteLine(title + " _sil_");
				for (int i = 0; i < lines.Length; i++) {
					silacWriters[silacIndex].WriteLine(lines[i]);
				}
				silacWriters[silacIndex].WriteLine("END IONS");
				silacWriters[silacIndex].WriteLine();
				silacCount++;
				return;
			}
			int isotopeIndex = peakList.GetIsotopeIndexForMsmsScanNumber(scanNumber);
			IsotopeCluster ic = null;
			if (isotopeIndex != -1) {
				ic = peakList.GetIsotopeCluster(isotopeIndex);
			}
			if (isotopeIndex != -1 && ic.IsotopeCorrelation >= isotopeCorrelationThreshold && ic.PolymerIndex == -1) {
				scanNumbers.Add(scanNumber);
				int ch = ic.Charge;
				double mass = peakList.GetMassEstimate(new int[] {isotopeIndex},
				                                       new int[] {ic.IsotopePatternStart}, new double[] {0}, false);
				double mz = mass / ch + MolUtil.MassProton;
				isoWriter.WriteLine("BEGIN IONS");
				isoWriter.WriteLine("PEPMASS=" + mz);
				isoWriter.WriteLine("CHARGE=" + ch + "+");
				isoWriter.WriteLine(title + " _iso_");
				for (int i = 0; i < lines.Length; i++) {
					isoWriter.WriteLine(lines[i]);
				}
				isoWriter.WriteLine("END IONS");
				isoWriter.WriteLine();
				isoCount++;
				return;
			}
			if (isotopeIndex != -1 && ic.PolymerIndex != -1) {
				scanNumbers.Add(scanNumber);
				int ch = ic.Charge;
				double mass = peakList.GetMassEstimate(new int[] {isotopeIndex},
				                                       new int[] {ic.IsotopePatternStart}, new double[] {0}, false);
				double mz = mass / ch + MolUtil.MassProton;
				polyWriter.WriteLine("BEGIN IONS");
				polyWriter.WriteLine("PEPMASS=" + mz);
				polyWriter.WriteLine("CHARGE=" + ch + "+");
				polyWriter.WriteLine(title + " _iso_");
				for (int i = 0; i < lines.Length; i++) {
					polyWriter.WriteLine(lines[i]);
				}
				polyWriter.WriteLine("END IONS");
				polyWriter.WriteLine();
				polyCount++;
				return;
			}
			scanNumbers.Add(scanNumber);
			peakWriter.WriteLine("BEGIN IONS");
			peakWriter.WriteLine(origMass);
			peakWriter.WriteLine("CHARGE=2+ and 3+");
			peakWriter.WriteLine(title + " _nix_");
			for (int i = 0; i < lines.Length; i++) {
				peakWriter.WriteLine(lines[i]);
			}
			peakWriter.WriteLine("END IONS");
			peakWriter.WriteLine();
			peakCount++;
		}

		private static int ExtractScanNumber(string title) {
			int index = title.IndexOf("FinneganScanNumber:");
			return int.Parse(title.Substring(index + "FinneganScanNumber:".Length).Trim());
		}

		public static int[] ExtractMsm(string filename, IPeakList peakList, IRawFile rawFile, string msmName,
		                               string msmNameBinary, string msmNameIndex, int topx,
		                               double isotopeCorrelationThreshold, SilacType type, out long[] pointers) {
			StreamWriter sw = new StreamWriter(msmName);
			BinaryWriter writer = FileUtil.GetBinaryWriter(msmNameBinary);
			List<long> filePointers = new List<long>();
			List<int> fileScanNumbers = new List<int>();
			long filePos = 0;
			for (int i = 0; i < rawFile.MS2Count; i++) {
				Spectrum s = rawFile.GetMS2Spectrum(i);
				int scanNumber = rawFile.GetScanNumberFromMs2Index(i);
				double simpleMz = rawFile.GetMS2MonoisotopicMz(i);
				if (double.IsNaN(simpleMz) || simpleMz < 1) {
					simpleMz = peakList.GetMs2Mz(i);
					if (double.IsNaN(simpleMz) || simpleMz < 1) {
						simpleMz = 1;
					}
				}
				int charge;
				double mz;
				double mass;
				GetPrecursorInfo(out charge, out mz, out mass, peakList, scanNumber, type, isotopeCorrelationThreshold, simpleMz);				
				if (rawFile.GetMS2SignalType(i) != SignalType.Centroid) {
					double[] specMasses;
					float[] specIntensities;
					T01PeakDetection.DetectPeaks(s, false, 3, CentroidPosition.gaussian, out specMasses, out specIntensities);
					s = new Spectrum(specMasses, specIntensities);
				}
				if (topx > 0) {
					s = s.TopX(topx, 100);
				}
				int sn = rawFile.GetScanNumberFromMs2Index(i);
				if (s.Count == 0) {
					continue;
				}
				StringBuilder sb = new StringBuilder();
				writer.Write(s.Count);
				for (int j = 0; j < s.Count; j++) {
					sb.AppendFormat("{0:F3}\t{1:F2}\n", s.GetMass(j), s.GetIntensity(j));
					writer.Write(s.GetMass(j));
					writer.Write(1.0 * s.GetIntensity(j));
				}
				filePointers.Add(filePos);
				fileScanNumbers.Add(sn);
				filePos += 16 * s.Count + 4;
				sw.WriteLine("BEGIN IONS");
				sw.WriteLine("PEPMASS={0:F6}", mz);
				if (charge > 0) {
					sw.WriteLine("CHARGE={0}+", charge);
				} else {
					sw.WriteLine("CHARGE=2+ and 3+");
				}
				sw.WriteLine("TITLE=RawFile: {0} FinneganScanNumber: {1}", filename, sn);
				sw.Write(sb.ToString());
				sw.WriteLine("END IONS");
				sw.WriteLine();
			}
			writer.Close();
			sw.Close();
			pointers = filePointers.ToArray();
			return fileScanNumbers.ToArray();
		}

		public static void GetPrecursorInfo(out int charge, out double mz, out double mass, IPeakList peakList, int scanNumber, 
			SilacType type, double isotopeCorrelationThreshold, double simpleMz) {
			charge = 0;
			mass = 0;
			int[] silacInfo = peakList.GetSilacInfoForMsmsScanNumber(scanNumber);
			if (silacInfo != null) {
				int silacId = silacInfo[0];
				int silacIndex = silacInfo[1];
				charge = peakList.GetSilacCharge(silacId);
				mass = peakList.GetSilacMass(silacId, silacIndex, type);
				mz = mass / charge + MolUtil.MassProton;
			} else {
				int isotopeIndex = peakList.GetIsotopeIndexForMsmsScanNumber(scanNumber);
				IsotopeCluster ic = null;
				if (isotopeIndex != -1) {
					ic = peakList.GetIsotopeCluster(isotopeIndex);
				}
				if (isotopeIndex != -1 && ic.IsotopeCorrelation >= isotopeCorrelationThreshold) {
					charge = ic.Charge;
					mass = peakList.GetMassEstimate(new int[] { isotopeIndex },
													new int[] { ic.IsotopePatternStart }, new double[] { 0 }, false);
					mz = mass / charge + MolUtil.MassProton;
				} else{
					int peakIndex = peakList.GetPeakIndexForMsmsScanNumber(scanNumber);
					mz = peakIndex >= 0 ? peakList.GetMz(peakIndex) : simpleMz;
				}
			}
		}

		public static void WriteSearchEngineParams(SearchEngineParams msh, string[] silacMsmNames, string isoMsmName,
		                                           string peakMsmName,
		                                           SilacLabel[][] l, SilacLabel[] allLabels, bool labelModsForNonSilac) {
			WriteSearchEngineParams(peakMsmName, allLabels, labelModsForNonSilac, true, (SearchEngineParams)msh.Clone(), " - " + ".peak.msm");
			WriteSearchEngineParams(isoMsmName, allLabels, labelModsForNonSilac, true, (SearchEngineParams)msh.Clone(), " - " + ".peak.msm");
			for (int i = 0; i < silacMsmNames.GetLength(0); i++) {
				WriteSearchEngineParams(silacMsmNames[i], l[i], true, false, (SearchEngineParams)msh.Clone(), " - " + ".sil" + i + ".msm");
			}
		}

		public static void WriteSearchEngineParams(string msmName, SilacLabel[] labels, bool addLabels, bool variableMods,
		                                           SearchEngineParams msh, string titleSuffix) {
			if (addLabels) {
				for (int i = 0; i < labels.Length; i++) {
					if (variableMods) {
						msh.AddVariabeModification(AminoAcid.GetMascotModificationStringForLabel(labels[i]));
					} else {
						msh.AddFixedModification(AminoAcid.GetMascotModificationStringForLabel(labels[i]));
					}
				}
			}
			msh.Title += titleSuffix;
			msh.Write(msmName.Substring(0, msmName.Length - 4) + ".par");
		}
	}
}