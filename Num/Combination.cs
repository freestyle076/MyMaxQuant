/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;

namespace MaxQuant.Num {
	internal class Combination {
		private readonly int n;
		private readonly int k;
		private readonly int[] data;

		internal Combination(int n, int k) {
			if (n < 0 || k < 0) {
				throw new Exception("Negative parameter in constructor");
			}
			this.n = n;
			this.k = k;
			data = new int[k];
			for (int i = 0; i < k; ++i) {
				data[i] = i;
			}
		}

		internal int[] Data {
			get { return data; }
		}

		internal Combination Successor {
			get {
				if (data[0] == n - k) {
					return null;
				}
				Combination ans = new Combination(n, k);
				long i;
				for (i = 0; i < k; ++i) {
					ans.data[i] = data[i];
				}
				for (i = k - 1; i > 0 && ans.data[i] == n - k + i; --i) {
				}
				++ans.data[i];
				for (long j = i; j < k - 1; ++j) {
					ans.data[j + 1] = ans.data[j] + 1;
				}
				return ans;
			}
		}
	}
}