/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Text;

namespace MaxQuant.Util{
	/// <summary>
	/// Static class containing utility routines for the manipulation of strings.
	/// </summary>
	public static class StringUtil{
		/// <summary>
		/// The digits 0 to 9 as subscripts.
		/// </summary>
		private static readonly char[] subscripts =
			new char[]{
			          	'\u2080', '\u2081', '\u2082', '\u2083', '\u2084', '\u2085', '\u2086', '\u2087', '\u2088',
			          	'\u2089'
			          };

		/// <summary>
		/// The digits 0 to 9 as superscripts.
		/// </summary>
		private static readonly char[] superscripts =
			new char[]{
			          	'\u2070', '\u00b9', '\u00b2', '\u00b3', '\u2074', '\u2075', '\u2076', '\u2077', '\u2078',
			          	'\u2079'
			          };

		/// <summary>
		/// Returns a string that is the same as the input string, except that all whitespace characters are removed.
		/// </summary>
		public static string RemoveWhitespace(string str){
			StringBuilder s = new StringBuilder();
			int len = str.Length;
			for (int i = 0; i < len; i++){
				char c = str[i];
				if (!Char.IsWhiteSpace(c)){
					s.Append(c);
				}
			}
			return s.ToString();
		}

		/// <summary>
		/// Returns a string that is the same as the input string, except that all consecutive sets of whitespace 
		/// characters are replaced by a single blank character.
		/// </summary>
		public static string ReduceWhitespace(string str){
			StringBuilder s = new StringBuilder();
			int len = str.Length;
			bool previousWasWhiteSpace = false;
			for (int i = 0; i < len; i++){
				char c = str[i];
				if (!Char.IsWhiteSpace(c)){
					s.Append(c);
					previousWasWhiteSpace = false;
				} else{
					if (!previousWasWhiteSpace){
						s.Append(' ');
					}
					previousWasWhiteSpace = true;
				}
			}
			return s.ToString().Trim();
		}

		/// <summary>
		/// Concatenates the string representations of the objects in the given array using the specified separator.
		/// </summary>
		/// <typeparam name="T">Type of objects to be concatenated as strings.</typeparam>
		/// <param name="separator">A string used to separate the array members.</param>
		/// <param name="o">The list of objects to be concatenated.</param>
		/// <returns>The concatenated string of all string representations of the array members.</returns>
		public static string Concat<T>(string separator, T[] o){
			if (o.Length == 0){
				return "";
			}
			if (o.Length == 1){
				return o[0].ToString();
			}
			StringBuilder s = new StringBuilder(o[0].ToString());
			for (int i = 1; i < o.Length; i++){
				s.Append(separator + o[i]);
			}
			return s.ToString();
		}

		/// <summary>
		/// Returns a string containing a representation of the given integer as superscript.
		/// </summary>
		/// <param name="n">The integer to be converted to superscript.</param>
		/// <param name="explicitPlus">Whether or not a '+' is added in front of positive numbers.</param>
		/// <returns>Representation of the given integer as superscript string.</returns>
		public static string ToSuperscript(int n, bool explicitPlus){
			bool isNegative = n < 0;
			bool isPositive = n > 0;
			n = Math.Abs(n);
			string nn = n.ToString();
			StringBuilder result = new StringBuilder();
			if (isNegative){
				result.Append('\u207B');
			}
			if (isPositive && explicitPlus){
				result.Append('\u207A');
			}
			char[] nnn = nn.ToCharArray();
			for (int i = 0; i < nnn.Length; i++){
				result.Append(superscripts[nnn[i] - '0']);
			}
			return result.ToString();
		}

		/// <summary>
		/// Returns a string containing a representation of the given integer as subscript.
		/// </summary>
		/// <param name="n">The integer to be converted to subscript.</param>
		/// <param name="explicitPlus">Whether or not a '+' is added in front of positive numbers.</param>
		/// <returns>Representation of the given integer as subscript string.</returns>
		public static string ToSubscript(int n, bool explicitPlus){
			bool isNegative = n < 0;
			bool isPositive = n > 0;
			n = Math.Abs(n);
			string nn = n.ToString();
			StringBuilder result = new StringBuilder();
			if (isNegative){
				result.Append('\u208B');
			}
			if (isPositive && explicitPlus){
				result.Append('\u208A');
			}
			char[] nnn = nn.ToCharArray();
			for (int i = 0; i < nnn.Length; i++){
				result.Append(subscripts[nnn[i] - '0']);
			}
			return result.ToString();
		}

		public static string Backslashify(string s, char c){
			StringBuilder result = new StringBuilder();
			for (int i = 0; i < s.Length; i++){
				char x = s[i];
				if (x != c){
					result.Append(x);
				} else{
					if (i == 0 || s[i - 1] != '\\'){
						result.Append("\\" + x);
					} else{
						result.Append(x);
					}
				}
			}
			return result.ToString();
		}

		public static int OccurenceCount(string s, char c){
			if(s == null){
				return 0;
			}
			int count = 0;
			foreach(char w in s){
				if(w == c){
					count++;
				}
			}
			return count;
		}
	}
}