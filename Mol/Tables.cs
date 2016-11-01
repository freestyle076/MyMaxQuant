/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Mol {
	public class Tables {
		public static readonly Dictionary<char, double> aaMonoMasses = InitMasses(); //TODO get rid of this

		public static readonly Dictionary<string, SequenceDatabase> databases = new Dictionary<string, SequenceDatabase>();
		public static readonly Dictionary<string, Enzyme> enzymes = new Dictionary<string, Enzyme>();
		public static readonly double MassNormalCTerminus = 15.994915 + 1.007825;
		public static readonly double MassNormalNTerminus = 1.007825;

		public static readonly Dictionary<string, Modification> modifications = new Dictionary<string, Modification>();
		public static Modification[] modificationList;

		/// <summary>
		/// Load the mascot configuration files.
		/// </summary>
		static Tables() {
			FileUtil.GetConfigPath();
			ParseDatabases();
			ParseEnzymes();
			ParseModifications();
		}

		private static Dictionary<char, double> InitMasses() {
			Dictionary<char, double> result = new Dictionary<char, double>();
			result.Add('A', 71.037114);
			result.Add('B', 114.534940);
			result.Add('C', 160.030649);
			result.Add('D', 115.026943);
			result.Add('E', 129.042593);
			result.Add('F', 147.068414);
			result.Add('G', 57.021464);
			result.Add('H', 137.058912);
			result.Add('I', 113.084064);
			result.Add('J', 0.000000);
			result.Add('K', 128.094963);
			result.Add('L', 113.084064);
			result.Add('M', 131.040485);
			result.Add('N', 114.042927);
			result.Add('O', 0.000000);
			result.Add('P', 97.052764);
			result.Add('Q', 128.058578);
			result.Add('R', 156.101111);
			result.Add('S', 87.032028);
			result.Add('T', 101.047679);
			result.Add('U', 150.953630);
			result.Add('V', 99.068414);
			result.Add('W', 186.079313);
			result.Add('X', 111.000000);
			result.Add('Y', 163.063329);
			result.Add('Z', 128.550590);
			return result;
		}

		public static string[] GetEnzymes() {
			string[] tmp = DataUtil.GetKeys(enzymes);
			Array.Sort(tmp);
			return tmp;
		}

		public static string[] GetDatabases() {
			string[] tmp = DataUtil.GetKeys(databases);
			Array.Sort(tmp);
			return tmp;
		}

		public static string[] GetModifications() {
			string[] result = DataUtil.GetKeys(modifications);
			Array.Sort(result);
			return result;
		}

		private static void ParseModifications() {
			StreamReader sr = new StreamReader(FileUtil.GetConfigPath() + "mod_file");
			string line;
			Modification modification = new Modification();
			List<char> aaRes = new List<char>();
			List<double> aaMass = new List<double>();
			List<double> aaWeight = new List<double>();
			List<ModificationSiteType> aaTermType = new List<ModificationSiteType>();
			List<double> neutralLoss = new List<double>();
			List<double> pepNeutralLoss = new List<double>();
			List<double> neutralLossWeight = new List<double>();
			List<double> pepNeutralLossWeight = new List<double>();
			List<Modification> mod = new List<Modification>();
			while ((line = sr.ReadLine()) != null) {
				if (line.Length > 0) {
					if (line[0] == '*') {
						modification.Site = aaRes.ToArray();
						if (aaMass.Count > 0) {
							modification.DeltaMass = aaMass[0];
							modification.SetPosition(ModificationPosition.anywhere);
						}
						modification.SetTermType(aaTermType.ToArray());
						modification.SetNeutralLoss(neutralLoss.ToArray());
						mod.Add(modification);
						modification = new Modification();
						aaRes = new List<char>();
						aaMass = new List<double>();
						aaWeight = new List<double>();
						aaTermType = new List<ModificationSiteType>();
						neutralLoss = new List<double>();
						pepNeutralLoss = new List<double>();
						neutralLossWeight = new List<double>();
						pepNeutralLossWeight = new List<double>();
					} else {
						int index = line.IndexOf(':');
						if (index < 0) {
							continue;
						}
						string[] m = new string[2];
						m[0] = line.Substring(0, index);
						m[1] = line.Substring(index + 1);
						m[1] = Regex.Replace(m[1], @"^\s+", "");
						string[] mm = Regex.Split(m[1], @"\s+");
						switch (m[0].ToLower()) {
							case "title":
								modification.Title = m[1];
								break;
							case "residues": {
								aaRes.Add(Char.ToUpper(Convert.ToChar(mm[0])));
								aaMass.Add(Convert.ToDouble(mm[1]));
								aaWeight.Add(Convert.ToDouble(mm[2]));
								aaTermType.Add(ModificationSiteType.aa);
								break;
							}
							case "residuesnterm": {
								aaRes.Add(Char.ToUpper(Convert.ToChar(mm[0])));
								aaMass.Add(Convert.ToDouble(mm[1]));
								aaWeight.Add(Convert.ToDouble(mm[2]));
								aaTermType.Add(ModificationSiteType.nterm);
								break;
							}
							case "residuescterm": {
								aaRes.Add(Char.ToUpper(Convert.ToChar(mm[0])));
								aaMass.Add(Convert.ToDouble(mm[1]));
								aaWeight.Add(Convert.ToDouble(mm[2]));
								aaTermType.Add(ModificationSiteType.cterm);
								break;
							}
							case "proteinnterm": {
								modification.DeltaMass = Convert.ToDouble(mm[0]);
								modification.SetPosition(ModificationPosition.proteinNterm);
								break;
							}
							case "proteincterm": {
								modification.DeltaMass = Convert.ToDouble(mm[0]);
								modification.SetPosition(ModificationPosition.proteinCterm);
								break;
							}
							case "neutralloss": {
								if (Convert.ToDouble(mm[0]) > 0.001) {
									neutralLoss.Add(Convert.ToDouble(mm[0]));
									neutralLossWeight.Add(Convert.ToDouble(mm[1]));
								}
								break;
							}
							case "pepneutralloss": {
								if (Convert.ToDouble(mm[0]) > 0.001) {
									pepNeutralLoss.Add(Convert.ToDouble(mm[0]));
									pepNeutralLossWeight.Add(Convert.ToDouble(mm[1]));
								}
								break;
							}
							case "nterm": {
								modification.DeltaMass = Convert.ToDouble(mm[0]);
								modification.SetPosition(ModificationPosition.anyNterm);
								break;
							}
							case "cterm": {
								modification.DeltaMass = Convert.ToDouble(mm[0]);
								modification.SetPosition(ModificationPosition.anyCterm);
								break;
							}
							default:
								break;
						}
					}
				}
			}
			sr.Close();
			Modification[] list = mod.ToArray();
			for (ushort i = 0; i < list.Length; i++) {
				if (list[i].AaCount > 0) {
					switch (list[i].GetTermType(0)) {
						case ModificationSiteType.aa:
							list[i].DeltaMass -= aaMonoMasses[list[i].GetAaAt(0)];
							break;
						case ModificationSiteType.nterm:
							list[i].DeltaMass -= MassNormalNTerminus;
							break;
						case ModificationSiteType.cterm:
							list[i].DeltaMass -= MassNormalCTerminus;
							break;
					}
				}
				if (list[i].GetPosition() == ModificationPosition.anyNterm) {
					list[i].DeltaMass -= MassNormalNTerminus;
				}
				if (list[i].GetPosition() == ModificationPosition.anyCterm) {
					list[i].DeltaMass -= MassNormalCTerminus;
				}
				if (list[i].GetPosition() == ModificationPosition.proteinNterm) {
					list[i].DeltaMass -= MassNormalNTerminus;
				}
				if (list[i].GetPosition() == ModificationPosition.proteinCterm) {
					list[i].DeltaMass -= MassNormalCTerminus;
				}
			}
			modificationList = mod.ToArray();
			for (ushort i = 0; i < modificationList.Length; i++) {
				modificationList[i].Index = i;
				modifications.Add(modificationList[i].Title, modificationList[i]);
			}
		}

		private static void ParseEnzymes() {
			StreamReader sr = new StreamReader(FileUtil.GetConfigPath() + "enzymes");
			string line;
			string title = null;
			Dictionary<int, string> cleavage = new Dictionary<int, string>();
			Dictionary<int, string> restrict = new Dictionary<int, string>();
			Dictionary<int, string> term = new Dictionary<int, string>();
			while ((line = sr.ReadLine()) != null) {
				if (line.Length > 0) {
					if (line[0] == '*') {
						int[] indices = DataUtil.GetKeys(cleavage);
						List<string> elements = new List<string>();
						foreach (int index in indices) {
							string cle = cleavage[index];
							string ter = term[index];
							string res = null;
							if (restrict.ContainsKey(index)) {
								res = restrict[index];
							}
							string elem;
							if (ter.Equals("C")) {
								if (res == null) {
									elem = "C." + cle;
								} else {
									elem = "C." + cle + "-" + res;
								}
							} else if (ter.Equals("N")) {
								if (res == null) {
									elem = "N." + cle;
								} else {
									elem = "N." + cle + "-" + res;
								}
							} else {
								throw new Exception("Cannot parse enzymes.");
							}
							elements.Add(elem);
						}
						enzymes.Add(title, new Enzyme(title, elements.ToArray()));
						title = null;
						cleavage = new Dictionary<int, string>();
						restrict = new Dictionary<int, string>();
						term = new Dictionary<int, string>();
					} else {
						string[] m = line.Split(':');
						string type = m[0];
						int index = 0;
						int x = type.IndexOf("[");
						int y = type.IndexOf("]");
						if (x != -1) {
							string s = type.Substring(x + 1, y - x - 1);
							index = int.Parse(s);
							type = type.Substring(0, x);
						}
						switch (type) {
							case "Title": {
								title = m[1];
								break;
							}
							case "Cleavage":
								cleavage.Add(index, m[1]);
								break;
							case "Restrict":
								restrict.Add(index, m[1]);
								break;
							case "Cterm":
								term.Add(index, "C");
								break;
							case "Nterm":
								term.Add(index, "N");
								break;
							default:
								break;
						}
					}
				}
			}
			sr.Close();
		}

		private static void ParseDatabases() {
			StreamReader sr = new StreamReader(FileUtil.GetConfigPath() + "mascot.dat");
			string line;
			List<string> buffer = null;
			string keyword = null;
			Dictionary<int, string> parseRules = new Dictionary<int, string>();
			List<string> databaseInfo = new List<string>();
			while ((line = sr.ReadLine()) != null) {
				line = line.Trim();
				if (line.StartsWith("#")) {
					continue;
				}
				if (line.ToLower().Equals("end")) {
					if (buffer == null) {
						return;
					}
					ParseBuffer(keyword, buffer.ToArray(), parseRules, databaseInfo);
					keyword = null;
					continue;
				}
				if (keyword == null) {
					keyword = line;
					buffer = new List<string>();
					continue;
				}
				buffer.Add(line);
			}
			sr.Close();
			foreach (string s in databaseInfo) {
				string[] w = s.Split('\t');
				string name = w[0];
				int parseRuleIndexAccession = Int32.Parse(w[1]);
				int parseRuleIndexDescription = Int32.Parse(w[2]);
				databases.Add(name,
				              new SequenceDatabase(name, parseRules[parseRuleIndexAccession], parseRules[parseRuleIndexDescription]));
			}
		}

		private static void ParseBuffer(string keyword, IEnumerable<string> lines, IDictionary<int, string> parseRules,
		                                ICollection<string> databaseInfo) {
			switch (keyword.ToLower().Trim()) {
				case "databases":
					foreach (string line in lines) {
						string[] w = line.Split('\t');
						databaseInfo.Add(w[0] + "\t" + w[10] + "\t" + w[11]);
					}
					break;
				case "parse":
					foreach (string line in lines) {
						string[] w = line.Split('\t');
						int id = Int32.Parse(w[0].Substring(5));
						string rule = w[1];
						if (rule.StartsWith("\"") && rule.EndsWith("\"")) {
							rule = rule.Substring(1, rule.Length - 2);
						}
						rule = rule.Replace("\\(", "(");
						rule = rule.Replace("\\)", ")");
						rule = StringUtil.Backslashify(rule, '|');
						parseRules.Add(id, rule);
					}
					break;
			}
		}

		public static Modification[] ToModifications(string[] modNames) {
			Modification[] result = new Modification[modNames.Length];
			for (int i = 0; i < result.Length; i++) {
				result[i] = modifications[modNames[i]];
			}
			return result;
		}

		public static ushort[] ToModificationIds(string[] modNames) {
			ushort[] result = new ushort[modNames.Length];
			for (int i = 0; i < result.Length; i++) {
				result[i] = modifications[modNames[i]].Index;
			}
			return result;
		}
	}
}