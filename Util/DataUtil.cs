/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Collections.Generic;

namespace MaxQuant.Util {

	public delegate void StringArg(string x);

	/// <summary>
	/// Static class containing utility routines for the manipulation of data structures.
	/// </summary>
	public static class DataUtil {
		/// <summary>
		/// Returns an array containing all keys from the given dictionary.
		/// </summary>
		public static K[] GetKeys<K, V>(Dictionary<K, V> dict) {
			K[] result = new K[dict.Count];
			dict.Keys.CopyTo(result, 0);
			return result;
		}

		/// <summary>
		/// Returns an array containing all values from the given dictionary.
		/// </summary>
		public static V[] GetValues<K, V>(Dictionary<K, V> dict) {
			V[] result = new V[dict.Count];
			dict.Values.CopyTo(result, 0);
			return result;
		}

		/// <summary>
		/// Returns an array containing all values from the given dictionary that correspond to a key in the given array.
		/// </summary>
		public static V[] GetValues<K, V>(Dictionary<K, V> dict, K[] keys) {
			V[] result = new V[keys.Length];
			for (int i = 0; i < result.Length; i++) {
				result[i] = dict[keys[i]];
			}
			return result;
		}
	}
}