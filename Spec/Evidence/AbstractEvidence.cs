/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.IO;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Spec.Evidence {
	public abstract class AbstractEvidence {
		private readonly int rawFileIndex;
		private List<MsmsHit> msmsHits = new List<MsmsHit>();
		private int id;
		private int modifiedPeptideId;
		private readonly int charge;
		private readonly float elutionTime;
		private readonly float calibratedElutionTime;

		protected AbstractEvidence(int rawFileIndex, float elutionTime, float calibratedElutionTime, int charge) {
			this.rawFileIndex = rawFileIndex;
			this.elutionTime = elutionTime;
			this.calibratedElutionTime = calibratedElutionTime;
			this.charge = charge;
		}

		protected AbstractEvidence(BinaryReader reader) {
			rawFileIndex = reader.ReadInt32();
			int len = reader.ReadInt32();
			for (int i = 0; i < len; i++) {
				msmsHits.Add(new MsmsHit(reader));
			}
			id = reader.ReadInt32();
			modifiedPeptideId = reader.ReadInt32();
			charge = reader.ReadInt32();
			elutionTime = reader.ReadSingle();
			calibratedElutionTime = reader.ReadSingle();
		}

		public int RawFileIndex {
			get { return rawFileIndex; }
		}

		public int Id {
			get { return id; }
			set { id = value; }
		}

		public int ModifiedPeptideId {
			get { return modifiedPeptideId; }
			set { modifiedPeptideId = value; }
		}

		public List<MsmsHit> MsmsHits {
			get { return msmsHits; }
		}

		public abstract string Type { get; }


		public float ElutionTime {
			get { return elutionTime; }
		}

		public float CalibratedElutionTime {
			get { return calibratedElutionTime; }
		}

		public int Charge {
			get { return charge; }
		}

		public int[] MsmsHitIds {
			get {
				List<MsmsHit> w = MsmsHits;
				int[] ids = new int[w.Count];
				for (int i = 0; i < ids.Length; i++) {
					ids[i] = w[i].Id;
				}
				return ids;
			}
		}

		public static AbstractEvidence Read(BinaryReader reader) {
			int which = reader.ReadInt32();
			switch (which) {
				case silacMsmsEvidence:
					return new SilacMsmsEvidence(reader);
				case silacPmfEvidence:
					return new SilacMassEvidence(reader);
				case SilacMatchEvidence:
					return new SilacMatchEvidence(reader);
				case isotopeMsmsEvidence:
					return new IsotopeMsmsEvidence(reader);
				case peakMsmsEvidence:
					return new PeakMsmsEvidence(reader);
			}
			throw new Exception("Never get here.");
		}

		public virtual void Write(BinaryWriter writer) {
			writer.Write(rawFileIndex);
			writer.Write(msmsHits.Count);
			for (int i = 0; i < msmsHits.Count; i++) {
				msmsHits[i].Write(writer);
			}
			writer.Write(id);
			writer.Write(modifiedPeptideId);
			writer.Write(charge);
			writer.Write(elutionTime);
			writer.Write(calibratedElutionTime);
		}

		public void AddMsmsHit(MsmsHit msmsHit) {
			msmsHits.Add(msmsHit);
		}

		public double GetPep() {
			double min = 1;
			for (int j = 0; j < msmsHits.Count; j++) {
				double fpr = msmsHits[j].Pep;
				if (fpr < min) {
					min = fpr;
				}
			}
			return min;
		}

		public double GetPep(int silacIndex) {
			double min = 1;
			for (int j = 0; j < msmsHits.Count; j++) {
				if (msmsHits[j].SilacIndex != silacIndex) {
					continue;
				}
				double fpr = msmsHits[j].Pep;
				if (fpr < min) {
					min = fpr;
				}
			}
			return min;
		}

		public MsmsHit GetBestMsmsHit(int silacIndex) {
			double min = double.MaxValue;
			int bestInd = -1;
			for (int j = 0; j < msmsHits.Count; j++) {
				if (msmsHits[j].SilacIndex != silacIndex) {
					continue;
				}
				double fpr = msmsHits[j].Pep;
				if (fpr < min) {
					min = fpr;
					bestInd = j;
				}
			}
			if (bestInd == -1) {
				return null;
			}
			return msmsHits[bestInd];
		}

		public double GetMascotScore() {
			MsmsHit h = GetBestMsmsHit();
			if (h == null) {
				return double.NaN;
			}
			return h.Score;
		}

		public double GetMascotScore(int silacIndex) {
			MsmsHit h = GetBestMsmsHit(silacIndex);
			if (h == null) {
				return 0;
			}
			return h.Score;
		}

		public double GetMascotDeltaScoreAll() {
			MsmsHit h = GetBestMsmsHit();
			if (h == null) {
				return double.NaN;
			}
			return h.DeltaScoreAll;
		}

		public double GetMascotDeltaScoreDifferentPep() {
			MsmsHit h = GetBestMsmsHit();
			if (h == null) {
				return double.NaN;
			}
			return h.DeltaScoreDifferentPep;
		}

		public MsmsHit GetBestMsmsHit() {
			int ind = -1;
			double max = -Double.MaxValue;
			for (int j = 0; j < msmsHits.Count; j++) {
				double score = msmsHits[j].PtmScore;
				if (score > max) {
					max = score;
					ind = j;
				}
			}
			if (ind == -1) {
				if (msmsHits.Count > 0) {
					return msmsHits[0];
				}
				return null;
			}
			return msmsHits[ind];
		}

		protected const int silacMsmsEvidence = 0;
		protected const int silacPmfEvidence = 1;
		protected const int SilacMatchEvidence = 2;
		protected const int isotopeMsmsEvidence = 3;
		protected const int peakMsmsEvidence = 4;

		public double GetPtmScore() {
			MsmsHit h = GetBestMsmsHit();
			if (h == null) {
				return double.NaN;
			}
			return h.PtmScore;
		}

		public double GetPtmDelta() {
			MsmsHit h = GetBestMsmsHit();
			if (h == null) {
				return double.NaN;
			}
			return h.PtmDelta;
		}

		public long GetPtmCombinatorics() {
			MsmsHit h = GetBestMsmsHit();
			if (h == null) {
				return 0;
			}
			return h.PtmCombinatorics;
		}

		public virtual double GetMassErrorPpm() {
			MsmsHit h = GetBestMsmsHit();
			if (h == null) {
				return double.NaN;
			}
			return h.MassErrorPpm;
		}

		public int[] GetSilacIndices() {
			HashSet<int> result = new HashSet<int>();
			foreach (MsmsHit msmsHit in msmsHits) {
				int ind = msmsHit.SilacIndex;
				if (ind >= 0 && !result.Contains(ind)) {
					result.Add(ind);
				}
			}
			int[] x = result.ToArray();
			Array.Sort(x);
			return x;
		}

		public virtual void Dispose() {
			if (msmsHits != null) {
				foreach (MsmsHit msmsHit in msmsHits) {
					msmsHit.Dispose();
				}
				msmsHits.Clear();
				msmsHits = null;
			}
		}
	}
}