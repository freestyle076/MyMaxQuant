/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public class SpectrumCache {
		private const int capacity = 1024 * 1024 * 250;

		private readonly Dictionary<int, Spectrum> cache = new Dictionary<int, Spectrum>();
		private int size;

		public Spectrum this[int index] {
			get { return cache[index]; }
		}

		public bool ContainsScanIndex(int index) {
			return cache.ContainsKey(index);
		}

		public void Add(int index, Spectrum spectrum) {
			if (size > capacity) {
				int[] nums = DataUtil.GetKeys(cache);
				Array.Sort(nums);
				nums = ArrayUtil.SubArray(nums, Math.Max(1, nums.Length / 10));
				foreach (int num in nums) {
					size -= 16 + cache[num].SizeEstimate;
					cache[num].Dispose();
					cache.Remove(num);
				}
			}
			cache.Add(index, spectrum);
			size += 16 + spectrum.SizeEstimate;
		}
	}
}