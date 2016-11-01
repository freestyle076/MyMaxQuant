/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.IO;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public class SearchEngineParams : ICloneable{
		private readonly string database;
		private readonly string enzyme;
		private List<string> fixedModifications = new List<string>();
		private readonly double fragmentIonTol;
		private readonly string fragmentIonTolUnit;
		private readonly string instrument;
		private readonly string massType;
		private readonly int maxMissedCleavages;
		private readonly double peptideMassTol;
		private readonly string peptideMassTolUnit;
		private readonly string searchType;
		private readonly string taxonomy;
		private readonly string username;
		private List<string> variableModifications = new List<string>();
		private string title;

		public SearchEngineParams(string enzyme, string database, double fragmentIonTol, string fragmentIonTolUnit,
		                          double peptideMassTol, string peptideMassTolUnit, string taxonomy, int maxMissedCleavages,
		                          IEnumerable<string> variableModifications, IEnumerable<string> fixedModifications,
		                          string massType, string searchType, string username, string title, string instrument) {
			this.enzyme = enzyme;
			this.database = database;
			this.fragmentIonTol = fragmentIonTol;
			this.fragmentIonTolUnit = fragmentIonTolUnit;
			this.peptideMassTol = peptideMassTol;
			this.peptideMassTolUnit = peptideMassTolUnit;
			this.taxonomy = taxonomy;
			this.maxMissedCleavages = maxMissedCleavages;
			foreach (string s in variableModifications) {
				this.variableModifications.Add(s);
			}
			foreach (string s in fixedModifications) {
				this.fixedModifications.Add(s);
			}
			this.massType = massType;
			this.searchType = searchType;
			this.username = username;
			this.title = title;
			this.instrument = instrument;
		}

		public SearchEngineParams(string filename) {
			StreamReader reader = new StreamReader(filename);
			string line;
			while ((line = reader.ReadLine()) != null) {
				line = line.Trim();
				if (line.Length == 0) {
					continue;
				}
				string[] w = line.Split('=');
				string key = w[0].Trim().ToUpper();
				string value = w[1].Trim();
				switch (key) {
					case "CLE":
						enzyme = value;
						break;
					case "DB":
						database = value;
						break;
					case "IT_MODS":
						variableModifications = new List<string>();
						if (value.Length > 0) {
							variableModifications.AddRange(value.Split(','));
						}
						break;
					case "ITOL":
						fragmentIonTol = double.Parse(value);
						break;
					case "ITOLU":
						fragmentIonTolUnit = value;
						break;
					case "MASS":
						massType = value;
						break;
					case "MODS":
						fixedModifications = new List<string>();
						if (value.Length > 0) {
							fixedModifications.AddRange(value.Split(','));
						}
						break;
					case "PFA":
						maxMissedCleavages = int.Parse(value);
						break;
					case "SEARCH":
						searchType = value;
						break;
					case "TAXONOMY":
						taxonomy = value;
						break;
					case "TITLE":
						title = value;
						break;
					case "TOL":
						peptideMassTol = double.Parse(value);
						break;
					case "TOLU":
						peptideMassTolUnit = value;
						break;
					case "USERNAME":
						username = value;
						break;
					case "INSTRUMENT":
						instrument = value;
						break;
					default:
						throw new Exception("Unknown key: " + key);
				}
			}
			reader.Close();
		}

		public string Title {
			get { return title; }
			set { title = value; }
		}

		public string Database {
			get { return database; }
		}

		public string[] FixedModifications {
			get { return fixedModifications.ToArray(); }
		}

		public string[] VariableModifications {
			get { return variableModifications.ToArray(); }
		}

		public int MaxMissedCleavages {
			get { return maxMissedCleavages; }
		}

		public string Enzyme {
			get { return enzyme; }
		}

		public void Write(string parameterFile) {
			StreamWriter sw = File.CreateText(parameterFile);
			sw.WriteLine("CLE={0}", enzyme);
			sw.WriteLine("DB={0}", database);
			sw.WriteLine("IT_MODS={0}", StringUtil.Concat(",", variableModifications.ToArray()));
			sw.WriteLine("ITOL={0}", fragmentIonTol);
			sw.WriteLine("ITOLU={0}", fragmentIonTolUnit);
			sw.WriteLine("MASS={0}", massType);
			sw.WriteLine("MODS={0}", StringUtil.Concat(",", fixedModifications.ToArray()));
			sw.WriteLine("PFA={0}", maxMissedCleavages);
			sw.WriteLine("SEARCH={0}", searchType);
			sw.WriteLine("TAXONOMY={0}", taxonomy);
			sw.WriteLine("TITLE={0}", title);
			sw.WriteLine("TOL={0}", peptideMassTol);
			sw.WriteLine("TOLU={0}", peptideMassTolUnit);
			sw.WriteLine("USERNAME={0}", username);
			sw.WriteLine("INSTRUMENT={0}", instrument);
			sw.Close();
		}

		public void AddVariabeModification(string s) {
			variableModifications.Add(s);
		}

		public void AddFixedModification(string s) {
			fixedModifications.Add(s);
		}

		public object Clone(){
			SearchEngineParams result = (SearchEngineParams)MemberwiseClone();
			result.fixedModifications = new List<string>(fixedModifications);
			result.variableModifications = new List<string>(variableModifications);
			return result;
		}
	}
}