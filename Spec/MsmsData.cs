/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.IO;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public class MsmsData {

		private string filename;
		private string indexName;
		private long[] filePointers;
		private int[] scanNumbers;
		private BinaryReader spectrumReader;

		public MsmsData(string filename, string indexName, long[] filePointers, int[] scanNumbers) : this(filename, indexName){
			this.filePointers = filePointers;
			this.scanNumbers = scanNumbers;
			WriteIndex();
		}

		public MsmsData(string filename, string indexName) {
			this.filename = filename;
			this.indexName = indexName;
		}

		public int[] ScanNumbers {
			get { return scanNumbers; }
		}

		public void WriteIndex() {
			BinaryWriter writer = FileUtil.GetBinaryWriter(indexName);
			writer.Write(filePointers.Length);
			for (int i = 0; i < filePointers.Length; i++) {
				writer.Write(filePointers[i]);
				writer.Write(scanNumbers[i]);
			}
			writer.Close();
		}

		public void ReadIndex() {
			if(!File.Exists(indexName)) {
				return;
			}
			BinaryReader reader = FileUtil.GetBinaryReader(indexName);
			int len = reader.ReadInt32();
			filePointers = new long[len];
			scanNumbers = new int[len];
			for (int i = 0; i < filePointers.Length; i++) {
				filePointers[i] = reader.ReadInt64();
				scanNumbers[i] = reader.ReadInt32();
			}
			reader.Close();
		}

		public void GetSpectrumFromScanNumber(int scanNumber, out double[] masses, out float[] intensities) {
			if (scanNumbers == null) {
				ReadIndex();
			}
			if (scanNumbers == null) {
				masses = null;
				intensities = null;
				return;
			}
			int ind = Array.BinarySearch(scanNumbers, scanNumber);
			if(ind < 0) {
				masses = null;
				intensities = null;
				return;
			}
			if (spectrumReader == null) {
				spectrumReader = FileUtil.GetBinaryReader(filename);
			}
			spectrumReader.BaseStream.Seek(filePointers[ind], SeekOrigin.Begin);
			int len = spectrumReader.ReadInt32();
			masses = new double[len];
			intensities = new float[len];
			for(int i = 0; i < len; i++) {
				masses[i] = spectrumReader.ReadDouble();
				intensities[i] = (float)spectrumReader.ReadDouble();
			}
		}

		public void Dispose() {
			if(spectrumReader != null) {
				spectrumReader.Close();
			}
			spectrumReader = null;
			filename = null;
			indexName = null;
			filePointers = null;
			scanNumbers = null;
		}
	}
}