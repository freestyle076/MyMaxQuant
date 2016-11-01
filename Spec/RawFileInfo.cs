/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Collections.Generic;
using System.IO;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public class RawFileInfo {
		private readonly string[] filePaths;
		private readonly string[] fileNames;
		private readonly Dictionary<string, ushort> gelSlices = new Dictionary<string, ushort>();
		private readonly HashSet<string> invert0 = new HashSet<string>();
		private readonly HashSet<string> invert1 = new HashSet<string>();
		private readonly HashSet<string> invert2 = new HashSet<string>();
		private readonly List<string> columnNames = new List<string>();
		private readonly List<Dictionary<string, string>> values = new List<Dictionary<string,string>>();

		public RawFileInfo(string[] filePaths) {
			this.filePaths = filePaths;
			fileNames = CreateFileNames(filePaths);
		}

		public RawFileInfo(BinaryReader reader) {
			int len = reader.ReadInt32();
			filePaths = new string[len];
			for (int i = 0; i < len; i++) {
				filePaths[i] = FileUtil.ReadString(reader);
			}
			len = reader.ReadInt32();
			fileNames = new string[len];
			for (int i = 0; i < len; i++) {
				fileNames[i] = FileUtil.ReadString(reader);
			}
			len = reader.ReadInt32();
			for (int i = 0; i < len; i++) {
				string key = FileUtil.ReadString(reader);
				ushort value = reader.ReadUInt16();
				gelSlices.Add(key, value);
			}
			len = reader.ReadInt32();
			for (int i = 0; i < len; i++) {
				string key = FileUtil.ReadString(reader);
				invert0.Add(key);
			}
			len = reader.ReadInt32();
			for (int i = 0; i < len; i++) {
				string key = FileUtil.ReadString(reader);
				invert1.Add(key);
			}
			len = reader.ReadInt32();
			for (int i = 0; i < len; i++) {
				string key = FileUtil.ReadString(reader);
				invert2.Add(key);
			}
			len = reader.ReadInt32();
			for (int i = 0; i < len; i++) {
				string name = FileUtil.ReadString(reader);
				columnNames.Add(name);
			}
			len = reader.ReadInt32();
			for (int i = 0; i < len; i++ ) {
				Dictionary<string, string> g = new Dictionary<string, string>();
				int len2 = reader.ReadInt32();
				for (int j = 0; j < len2; j++) {
					string key = FileUtil.ReadString(reader);
					string value = FileUtil.ReadString(reader);
					g.Add(key, value);
				}
				values.Add(g);
			}
		}

		public int Count {
			get { return filePaths.Length; }
		}

		public bool HasGelSlices {
			get { return gelSlices.Count > 0; }
		}

		public bool HasExperiment {
			get { return (GetExperimentIndex() != -1 && GetAllExperimentValues().Length > 0); }
		}

		private const string experimentName = "Experiment";

		public int GetExperimentIndex() {
			for (int i = 0; i < columnNames.Count; i++) {
				if (columnNames[i].Equals(experimentName)) {
					return i;
				}
			}
			return -1;
		}

		public string[] GetAllExperimentValues() {
			int index = GetExperimentIndex();
			Dictionary<string, string> v = values[index];
			string[] vals = DataUtil.GetValues(v);
			return ArrayUtil.UniqueValues(vals);
		}

		public void Write(BinaryWriter writer) {
			writer.Write(filePaths.Length);
			for (int i = 0; i < filePaths.Length; i++) {
				FileUtil.WriteString(filePaths[i], writer);
			}
			writer.Write(fileNames.Length);
			for (int i = 0; i < fileNames.Length; i++) {
				FileUtil.WriteString(fileNames[i], writer);
			}
			writer.Write(gelSlices.Count);
			string[] keys = DataUtil.GetKeys(gelSlices);
			foreach (string key in keys) {
				FileUtil.WriteString(key, writer);
				writer.Write(gelSlices[key]);
			}
			writer.Write(invert0.Count);
			foreach (string key in invert0.ToArray()) {
				FileUtil.WriteString(key, writer);
			}
			writer.Write(invert1.Count);
			foreach (string key in invert1.ToArray()) {
				FileUtil.WriteString(key, writer);
			}
			writer.Write(invert2.Count);
			foreach (string key in invert2.ToArray()) {
				FileUtil.WriteString(key, writer);
			}
			writer.Write(columnNames.Count);
			foreach(string name in columnNames) {
				FileUtil.WriteString(name, writer);
			}
			writer.Write(values.Count);
			foreach (Dictionary<string, string> v in values) {
				writer.Write(v.Count);
				keys = DataUtil.GetKeys(v);
				foreach (string key in keys) {
					FileUtil.WriteString(key, writer);
					FileUtil.WriteString(v[key], writer);
				}
			}
		}

		private static string[] CreateFileNames(string[] paths) {
			string[] result = new string[paths.Length];
			for (int i = 0; i < paths.Length; i++) {
				result[i] = CreateFileName(paths[i]);
			}
			return result;
		}

		private static string CreateFileName(string path) {
			int i0 = path.LastIndexOf("\\");
			int i1 = path.LastIndexOf(".");
			return path.Substring(i0 + 1, i1 - i0 - 1);
		}

		public void SetData(string fileName, string column, string value) {
			int index = columnNames.IndexOf(column);
			if(index == -1) {
				columnNames.Add(column);
				values.Add(new Dictionary<string, string>());
				index = columnNames.Count - 1;
			}
			values[index].Add(fileName, value);
		}

		public string GetFileName(int index) {
			return fileNames[index];
		}

		public ushort GetGelSlice(int index) {
			string fileName = fileNames[index];
			if(!gelSlices.ContainsKey(fileName)) {
				return ushort.MaxValue;
			}
			return gelSlices[fileName];
		}

		public int GetRawFileIndex(string name) {//TODO: use hash
			for(int i = 0; i < fileNames.Length; i++) {
				if(fileNames[i].Equals(name)) {
					return i;
				}
			}
			return -1;
		}

		public void SetSlice(string name, ushort slice) {
			gelSlices[name] = slice;
		}

		public bool HasFileName(string name) {
			foreach(string x in fileNames) {
				if(x.Equals(name)) {
					return true;
				}
			}
			return false;
		}

		public string[] GetColumnNames() {
			return columnNames.ToArray();
		}

		public ushort GetGelSlice(string name) {
			if(!gelSlices.ContainsKey(name)) {
				return ushort.MaxValue;
			}
			return gelSlices[name];
		}

		public string GetValue(int index, string colName) {
			return GetValue(fileNames[index], colName);
		}

		public string GetValue(string name, string colName) {
			int index = -1;
			for(int i  = 0; i < columnNames.Count; i++) {
				if(columnNames[i].Equals(colName)) {
					index = i;
					break;
				}
			}
			Dictionary<string, string> x = values[index];
			if(!x.ContainsKey(name)) {
				return null;
			}
			return x[name];
		}

		public ushort[] GetAllSlices() {
			ushort[] result = DataUtil.GetValues(gelSlices);
			return ArrayUtil.UniqueValues(result);
		}

		public string[] GetAllValues(string colName) {
			int index = columnNames.IndexOf(colName);
			Dictionary<string, string> v = values[index];
			string[] vals = DataUtil.GetValues(v);
			return ArrayUtil.UniqueValues(vals);
		}

		public HashSet<int> GetRawFileIndicesFromExperiment(string exp) {
			int index = GetExperimentIndex();
			Dictionary<string, string> v = values[index];
			HashSet<int> result = new HashSet<int>();
			for(int i = 0; i < fileNames.Length; i++) {
				string name = fileNames[i];
				if(v.ContainsKey(name)) {
					string val = v[name];
					if(val.Equals(exp)) {
						result.Add(i);
					}
				}
			}
			return result;
		}

		public string GetExperimentName(string path) {
			int index = GetExperimentIndex();
			Dictionary<string, string> v = values[index];
			string file = CreateFileName(path);
			if(!v.ContainsKey(file)) {
				return null;
			}
			return v[file];
		}

		public void SetInvert0(string name, bool invert) {
			if(invert) {
				invert0.Add(name);
			}
		}

		public void SetInvert1(string name, bool invert) {
			if (invert) {
				invert1.Add(name);
			}
		}

		public void SetInvert2(string name, bool invert) {
			if (invert) {
				invert2.Add(name);
			}
		}

		public string GetExperimentValue(int i) {
			return values[GetExperimentIndex()][fileNames[i]];
		}

		public bool Invert0(int index) {
			return invert0.Contains(fileNames[index]);
		}

		public bool Invert1(int index) {
			return invert1.Contains(fileNames[index]);
		}

		public bool Invert2(int index) {
			return invert2.Contains(fileNames[index]);
		}
	}
}