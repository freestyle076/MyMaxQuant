/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Collections.Generic;
using MaxQuant.Mol;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public class MascotQuery {
		private readonly Dictionary<int, MascotPeptide> peptides = new Dictionary<int, MascotPeptide>();

		private readonly int scanNum;

		private string rawFile;

		private readonly MascotQueryType queryType;

		private ushort[] fixedMods;

		public MascotQuery(int scanNum, string rawFile, MascotQueryType queryType, string[] fixedMods) {
			this.scanNum = scanNum;
			this.rawFile = rawFile;
			this.queryType = queryType;
			this.fixedMods = new ushort[fixedMods.Length];
			for(int i = 0; i < fixedMods.Length; i++) {
				this.fixedMods[i] = Tables.modifications[fixedMods[i]].Index;		
			}		
		}

		public void AddPeptide(int peptideNumber, MascotPeptide peptide) {
			peptides.Add(peptideNumber, peptide);
		}

		public MascotPeptide GetPeptide(int num) {
			if (!peptides.ContainsKey(num)) {
				return null;
			}
			return peptides[num];
		}

		public MascotPeptide[] GetPeptides() {
			return DataUtil.GetValues(peptides);
		}

		public string RawFile {
			get { return rawFile; }
			set { rawFile = value; }
		}

		public MascotQueryType QueryType {
			get { return queryType; }
		}

		public int ScanNum {
			get { return scanNum; }
		}

		public ushort[] FixedMods {
			get { return fixedMods; }
		}

		public void Dispose() {
			fixedMods = null;
			foreach (KeyValuePair<int, MascotPeptide> pair in peptides) {
				pair.Value.Dispose();
			}
			peptides.Clear();
		}

		public void RecalcMassDeviations(IPeakList peakList, int scanNumber, SilacType type, double isotopeCorrelationThreshold) {
			foreach (MascotPeptide mp in peptides.Values) {
				mp.RecalcMassDeviation(fixedMods, peakList, scanNumber, type, isotopeCorrelationThreshold);
			}
		}
	}
}