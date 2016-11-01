/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Collections;
using System.Collections.Generic;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public class NeighbourList {
		private readonly Dictionary<int, int[]> neighborListMultiple;
		private readonly Dictionary<int, int> neighborListSingle;

		public NeighbourList() {
			neighborListSingle = new Dictionary<int, int>();
			neighborListMultiple = new Dictionary<int, int[]>();
		}

		public void Add(int i, int j) {
			Add2(i, j);
			Add2(j, i);
		}

		private void Add2(int i, int j) {
			int d = j - i;

			if (neighborListMultiple.ContainsKey(i)) 
			{
				int[] w = neighborListMultiple[i];
				int[] w2 = ArrayUtil.Concat(w, d);
				neighborListMultiple.Remove(i);
				neighborListMultiple.Add(i, w2);
				return;
			}

			if (neighborListSingle.ContainsKey(i)) 
			{
				int w = neighborListSingle[i];
				int[] w2 = new int[] {w, d};
				neighborListSingle.Remove(i);
				neighborListMultiple.Add(i, w2);
				return;
			}

			neighborListSingle.Add(i, d);
		}

		public bool IsEmptyAt(int i) {
			if (neighborListSingle.ContainsKey(i)) {
				return false;
			}
			if (neighborListMultiple.ContainsKey(i)) {
				return false;
			}
			return true;
		}

		public IEnumerable GetAt(int i) {
			if (neighborListSingle.ContainsKey(i)) {
				return new int[] {neighborListSingle[i] + i};
			}
			if (neighborListMultiple.ContainsKey(i)) {
				int[] l = neighborListMultiple[i];
				int[] result = new int[l.Length];
				for (int j = 0; j < result.Length; j++) {
					result[j] = l[j] + i;
				}
				return result;
			}
			return null;
		}

		public void DeleteAt(int i) {
			if (neighborListSingle.ContainsKey(i)) {
				neighborListSingle.Remove(i);
				return;
			}
			if (neighborListMultiple.ContainsKey(i)) {
				neighborListMultiple.Remove(i);
			}
		}
	}
}