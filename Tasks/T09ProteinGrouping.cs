/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using MaxQuant.Util;

namespace MaxQuant.Tasks {
	public static class T09ProteinGrouping {
		public static void ClusterProteins(ref string[][] proteinIds, ref string[][] pepSeqs) {
			int n = proteinIds.Length;
			for (int i = 0; i < n; i++) {
				Array.Sort(pepSeqs[i]);
			}
			IndexedBitMatrix contains = new IndexedBitMatrix(n, n);
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i == j) {
						continue;
					}
					contains.Set(i, j, Contains(pepSeqs[i], pepSeqs[j]));
				}
			}
			int count;
			do {
				count = 0;
				int start = 0;
				while (true) {
					int container = -1;
					int contained = -1;
					for (int i = start; i < proteinIds.Length; i++) {
						container = GetContainer(i, contains);
						if (container != -1) {
							contained = i;
							break;
						}
					}
					if (container == -1) {
						break;
					}
					for (int i = 0; i < n; i++) {
						contains.Set(i, contained, false);
						contains.Set(contained, i, false);
					}
					pepSeqs[contained] = new string[0];
					proteinIds[container] = ArrayUtil.Concat(proteinIds[container], proteinIds[contained]);
					proteinIds[contained] = new string[0];
					start = contained + 1;
					count++;
				}
			} while (count > 0);
			List<int> valids = new List<int>();
			for (int i = 0; i < n; i++) {
				if (pepSeqs[i].Length > 0) {
					valids.Add(i);
				}
			}
			int[] a = valids.ToArray();
			proteinIds = ArrayUtil.SubArray(proteinIds, a);
			pepSeqs = ArrayUtil.SubArray(pepSeqs, a);
		}

		private static int GetContainer(int contained, IndexedBitMatrix contains) {
			int n = contains.RowCount;
			for (int i = 0; i < n; i++) {
				if (contains.Get(i, contained)) {
					return i;
				}
			}
			return -1;
		}

		private static bool Contains(string[] p1, ICollection<string> p2) {
			if (p2.Count > p1.Length) {
				return false;
			}
			foreach (string p in p2) {
				int index = Array.BinarySearch(p1, p);
				if (index < 0) {
					return false;
				}
			}
			return true;
		}
	}
}
