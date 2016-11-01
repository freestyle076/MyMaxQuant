/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.IO;
using System.IO.Compression;
using System.Reflection;
using System.Text;
using System.Windows.Forms;

namespace MaxQuant.Util {
	/// <summary>
	/// Static class containing utility routines for accessing and handling files.
	/// </summary>
	public class FileUtil {
		/// <summary>
		/// Creates a <code>BinaryReader</code> reading from the given file path.
		/// </summary>
		/// <param name="path">File to read from.</param>
		/// <returns>The <code>BinaryReader</code>.</returns>
		public static BinaryReader GetBinaryReader(string path) {
			return new BinaryReader(new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read));
		}

		public static StreamReader GetGzipReader(string filename) {
			return new StreamReader(new GZipStream(new FileStream(filename, FileMode.Open, FileAccess.Read), CompressionMode.Decompress));
		}

		/// <summary>
		/// Creates a <code>StreamReader</code> reading from the given text resource within this assembly.
		/// </summary>
		/// <param name="name">Name of the resource to read from.</param>
		/// <returns>The <code>StreamReader</code>.</returns>
		public static StreamReader GetResourceTextReader(string name) {
			return new StreamReader(GetResourceStream(name));
		}

		/// <summary>
		/// Creates a <code>Stream</code> reading from the given text resource within this assembly.
		/// </summary>
		/// <param name="name">Name of the resource to read from.</param>
		/// <returns>The <code>StreamReader</code>.</returns>
		public static Stream GetResourceStream(string name) {
			Assembly assembly = Assembly.GetExecutingAssembly();
			return assembly.GetManifestResourceStream(name);
		}

		/// <summary>
		/// Reads a string in a binary version which is purely ascii-encoded.
		/// </summary>
		public static string ReadString(BinaryReader reader) {
			int len = reader.ReadInt32();
			byte[] x = new byte[len];
			for (int i = 0; i < len; i++) {
				x[i] = reader.ReadByte();
			}
			char[] chars = new char[Encoding.ASCII.GetCharCount(x, 0, x.Length)];
			Encoding.ASCII.GetChars(x, 0, x.Length, chars, 0);
			return new string(chars);
		}

		/// <summary>
		/// Writes a string in a binary version which is purely ascii-encoded.
		/// </summary>
		public static void WriteString(string str, BinaryWriter writer) {
			byte[] x = Encoding.ASCII.GetBytes(str);
			writer.Write(x.Length);
			for (int i = 0; i < x.Length; i++) {
				writer.Write(x[i]);
			}
		}

		/// <summary>
		/// Deletes file after checking for its existence.
		/// </summary>
		public static void DeleteFile(string filename) {
			if (File.Exists(filename)) {
				File.Delete(filename);
			}
		}

		public static void WriteDoubleArray(double[] x, BinaryWriter writer) {
			writer.Write(x.Length);
			for (int i = 0; i < x.Length; i++) {
				writer.Write(x[i]);
			}
		}

		public static double[] ReadDoubleArray(BinaryReader reader) {
			int n = reader.ReadInt32();
			double[] result = new double[n];
			for (int i = 0; i < n; i++) {
				result[i] = reader.ReadDouble();
			}
			return result;
		}
		/// <summary>
		/// Removes all files and folders recursively in the specified folder 
		/// and the specified folder itself.
		/// </summary>
		/// <param name="path">Path of the folder to be removed.</param>
		public static void Rmdir(string path) {
			if (!Directory.Exists(path)) {
				return;
			}
			string[] d = Directory.GetDirectories(path);
			foreach (string s in d) {
				Rmdir(s);
			}
			string[] f = Directory.GetFiles(path);
			foreach (string s in f) {
				File.Delete(s);
			}
			Directory.Delete(path);
		}

		/// <summary>
		/// Creates a <code>BinaryWriter</code> writing to the given file path.
		/// </summary>
		/// <param name="path">File to write to.</param>
		/// <returns>The <code>BinaryWriter</code>.</returns>
		public static BinaryWriter GetBinaryWriter(string path) {
			return new BinaryWriter(new FileStream(path, FileMode.Create, FileAccess.Write));
		}

		public static string GetConfigPath() {
			string p = Application.ExecutablePath;
			p = p.Substring(0, p.LastIndexOf("\\")) + "\\conf\\";
			return p;
		}

		/// <summary>
		/// Tests whether the directory corresponding to the given path is writable.
		/// </summary>
		/// <param name="path">Path of the directory.</param>
		/// <returns><code>true</code> if the directory corresponding to the given path is writable.</returns>
		public static bool TestDirWritable(string path) {
			if (Directory.Exists(path)) {
				try {
					string FilePath = Path.Combine(path, "AccesTest.tmp");
					StreamWriter AccessTest = new StreamWriter(FilePath);
					AccessTest.Close();
					File.Delete(FilePath);
				} catch {
					return false;
				}
				return true;
			}
			return false;
		}
	}
}