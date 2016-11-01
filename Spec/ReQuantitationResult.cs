/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.IO;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public class ReQuantitationResult {
		private readonly string filename;
		private Dictionary<int, double[]> doubletRatioInfo;
		private Dictionary<int, double[][]> tripletRatioInfo;
		private Dictionary<int, double[]> doubletIntensityInfo;
		private Dictionary<int, double[]> tripletIntensityInfo;

		public ReQuantitationResult(string filename) {
			this.filename = filename;
		}

		public ReQuantitationResult(string filename, Dictionary<int, double[]> doubletRatioInfo,
		                            Dictionary<int, double[][]> tripletRatioInfo,
		                            Dictionary<int, double[]> doubletIntensityInfo,
		                            Dictionary<int, double[]> tripletIntensityInfo) {
			this.filename = filename;
			this.doubletRatioInfo = doubletRatioInfo;
			this.tripletRatioInfo = tripletRatioInfo;
			this.doubletIntensityInfo = doubletIntensityInfo;
			this.tripletIntensityInfo = tripletIntensityInfo;
			Write();
		}

		public double GetDoubletRatio(int index) {
			if (doubletRatioInfo == null) {
				Read();
			}
			if (doubletRatioInfo == null) {
				return double.NaN;
			}
			if (!doubletRatioInfo.ContainsKey(index)) {
				return double.NaN;
			}
			return doubletRatioInfo[index][0];
		}

		public double GetDoubletNormalizedRatio(int index) {
			if (doubletRatioInfo == null) {
				Read();
			}
			if (doubletRatioInfo == null) {
				return double.NaN;
			}
			if (!doubletRatioInfo.ContainsKey(index)) {
				return double.NaN;
			}
			return doubletRatioInfo[index][1];
		}

		public double GetDoubletRatioSignificanceA(int index) {
			if (doubletRatioInfo == null) {
				Read();
			}
			if (doubletRatioInfo == null) {
				return double.NaN;
			}
			if (!doubletRatioInfo.ContainsKey(index)) {
				return double.NaN;
			}
			return doubletRatioInfo[index][2];
		}

		public double GetDoubletRatioSignificanceB(int index) {
			if (doubletRatioInfo == null) {
				Read();
			}
			if (doubletRatioInfo == null) {
				return double.NaN;
			}
			if (!doubletRatioInfo.ContainsKey(index)) {
				return double.NaN;
			}
			return doubletRatioInfo[index][3];
		}

		public double[] GetDoubletIntensities(int index) {
			if (doubletRatioInfo == null) {
				Read();
			}
			if (doubletRatioInfo == null) {
				return null;
			}
			if (!doubletIntensityInfo.ContainsKey(index)) {
				return null;
			}
			return doubletIntensityInfo[index];
		}

		public double[] GetTripletIntensities(int index) {
			if (doubletRatioInfo == null) {
				Read();
			}
			if (doubletRatioInfo == null) {
				return null;
			}
			if (!tripletIntensityInfo.ContainsKey(index)) {
				return null;
			}
			return tripletIntensityInfo[index];
		}

		public double[] GetTripletRatio(int index) {
			if (tripletRatioInfo == null) {
				Read();
			}
			if (tripletRatioInfo == null) {
				return null;
			}
			if (!tripletRatioInfo.ContainsKey(index)) {
				return null;
			}
			return tripletRatioInfo[index][0];
		}

		public double[] GetTripletNormalizedRatio(int index) {
			if (tripletRatioInfo == null) {
				Read();
			}
			if (tripletRatioInfo == null) {
				return null;
			}
			if (!tripletRatioInfo.ContainsKey(index)) {
				return null;
			}
			return tripletRatioInfo[index][1];
		}

		public double[] GetTripletRatioSignificanceA(int index) {
			if (tripletRatioInfo == null) {
				Read();
			}
			if (tripletRatioInfo == null) {
				return null;
			}
			if (!tripletRatioInfo.ContainsKey(index)) {
				return null;
			}
			return tripletRatioInfo[index][2];
		}

		public double[] GetTripletRatioSignificanceB(int index) {
			if (tripletRatioInfo == null) {
				Read();
			}
			if (tripletRatioInfo == null) {
				return null;
			}
			if (!tripletRatioInfo.ContainsKey(index)) {
				return null;
			}
			return tripletRatioInfo[index][3];
		}

		public void Write() {
			BinaryWriter writer = FileUtil.GetBinaryWriter(filename);
			writer.Write(doubletRatioInfo.Count);
			foreach (int key in doubletRatioInfo.Keys) {
				writer.Write(key);
				double[] x = doubletRatioInfo[key];
				for (int i = 0; i < 4; i++) {
					writer.Write(x[i]);
				}
			}
			writer.Write(tripletRatioInfo.Count);
			foreach (int key in tripletRatioInfo.Keys) {
				writer.Write(key);
				double[][] x = tripletRatioInfo[key];
				for (int j = 0; j < 4; j++) {
					for (int k = 0; k < 3; k++) {
						writer.Write(x[j][k]);
					}
				}
			}
			writer.Write(doubletIntensityInfo.Count);
			foreach (int key in doubletIntensityInfo.Keys) {
				writer.Write(key);
				double[] x = doubletIntensityInfo[key];
				for (int i = 0; i < 2; i++) {
					writer.Write(x[i]);
				}
			}
			writer.Write(tripletIntensityInfo.Count);
			foreach (int key in tripletIntensityInfo.Keys) {
				writer.Write(key);
				double[] x = tripletIntensityInfo[key];
				for (int i = 0; i < 3; i++) {
					writer.Write(x[i]);
				}
			}
			writer.Close();
		}

		public void Read() {
			if (!File.Exists(filename)) {
				return;
			}
			doubletRatioInfo = new Dictionary<int, double[]>();
			tripletRatioInfo = new Dictionary<int, double[][]>();
			doubletIntensityInfo = new Dictionary<int, double[]>();
			tripletIntensityInfo = new Dictionary<int, double[]>();
			BinaryReader reader = FileUtil.GetBinaryReader(filename);
			int len = reader.ReadInt32();
			for (int i = 0; i < len; i++) {
				int key = reader.ReadInt32();
				double[] value = new double[4];
				for (int j = 0; j < 4; j++) {
					value[j] = reader.ReadDouble();
				}
				doubletRatioInfo.Add(key, value);
			}
			len = reader.ReadInt32();
			for (int i = 0; i < len; i++) {
				int key = reader.ReadInt32();
				double[][] value = new double[4][];
				for (int j = 0; j < 4; j++) {
					value[j] = new double[3];
					for (int k = 0; k < 3; k++) {
						value[j][k] = reader.ReadDouble();
					}
				}
				tripletRatioInfo.Add(key, value);
			}
			len = reader.ReadInt32();
			for (int i = 0; i < len; i++) {
				int key = reader.ReadInt32();
				double[] value = new double[2];
				for (int j = 0; j < 2; j++) {
					value[j] = reader.ReadDouble();
				}
				doubletIntensityInfo.Add(key, value);
			}
			len = reader.ReadInt32();
			for (int i = 0; i < len; i++) {
				int key = reader.ReadInt32();
				double[] value = new double[3];
				for (int j = 0; j < 3; j++) {
					value[j] = reader.ReadDouble();
				}
				tripletIntensityInfo.Add(key, value);
			}
			reader.Close();
		}

		public bool ContainsIsotopeIndex(int id, SilacType type) {
			if (doubletRatioInfo == null) {
				Read();
			}
			switch (type) {
				case SilacType.Doublets:
					return doubletRatioInfo.ContainsKey(id);
				case SilacType.Triplets:
					return tripletRatioInfo.ContainsKey(id);
				default:
					throw new Exception("Impossible.");
			}
		}

		public void Dispose() {
			if (doubletRatioInfo != null) {
				doubletRatioInfo.Clear();
				doubletRatioInfo = null;
			}
			if (tripletRatioInfo != null) {
				tripletRatioInfo.Clear();
				tripletRatioInfo = null;
			}
			if (doubletIntensityInfo != null) {
				doubletIntensityInfo.Clear();
				doubletIntensityInfo = null;
			}
			if (tripletIntensityInfo != null) {
				tripletIntensityInfo.Clear();
				tripletIntensityInfo = null;
			}
		}

		public void Add(Dictionary<int, double[]> doubletRatioInfoIn, Dictionary<int, double[][]> tripletRatioInfoIn,
		                Dictionary<int, double[]> doubletIntensityInfoIn, Dictionary<int, double[]> tripletIntensityInfoIn) {
			if (doubletRatioInfo == null) {
				Read();
			}
			foreach (KeyValuePair<int, double[]> pair in doubletRatioInfoIn) {
				if (!doubletRatioInfo.ContainsKey(pair.Key)) {
					doubletRatioInfo.Add(pair.Key, pair.Value);
				}
			}
			foreach (KeyValuePair<int, double[][]> pair in tripletRatioInfoIn) {
				if (!tripletRatioInfo.ContainsKey(pair.Key)) {
					tripletRatioInfo.Add(pair.Key, pair.Value);
				}
			}
			foreach (KeyValuePair<int, double[]> pair in doubletIntensityInfoIn) {
				if (!doubletIntensityInfo.ContainsKey(pair.Key)) {
					doubletIntensityInfo.Add(pair.Key, pair.Value);
				}
			}
			foreach (KeyValuePair<int, double[]> pair in tripletIntensityInfoIn) {
				if (!tripletIntensityInfo.ContainsKey(pair.Key)) {
					tripletIntensityInfo.Add(pair.Key, pair.Value);
				}
			}
		}
	}
}