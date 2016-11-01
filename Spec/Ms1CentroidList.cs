/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public class Ms1CentroidList {
		private double[] peakCenterMass;
		private float[] peakIntensity;
		private float[] peakMaxMass;
		private float[] peakMinMass;
		private GrowablePeak[] peaks;

		private Ms1CentroidList() {
		}

		public Ms1CentroidList(double[] peakCenterMass, float[] peakMinMass,
		                       float[] peakMaxMass, float[] peakIntensity) {
			this.peakCenterMass = peakCenterMass;
			this.peakMinMass = peakMinMass;
			this.peakMaxMass = peakMaxMass;
			this.peakIntensity = peakIntensity;
			peaks = new GrowablePeak[peakCenterMass.Length];
		}

		public int Count {
			get { return peakCenterMass.Length; }
		}

		public double CenterMass(int index) {
			if (index < 0) {
				return Double.NaN;
			}
			return peakCenterMass[index];
		}

		public float MaxMass(int index) {
			return peakMaxMass[index];
		}

		public GrowablePeak GetPeak(int index) {
			return peaks[index];
		}

		public float MinMass(int index) {
			return peakMinMass[index];
		}

		public float GetIntensity(int index) {
			return peakIntensity[index];
		}

		public int GetClosestIndex(double mass) {
			if (Count == 0) {
				return -1;
			}
			if (mass <= peakCenterMass[0]) {
				return 0;
			}
			if (mass >= peakCenterMass[Count - 1]) {
				return Count - 1;
			}
			int index = Array.BinarySearch(peakCenterMass, mass);
			if (index >= 0) {
				return index;
			}
			index = -2 - index;
			if (Math.Abs(peakCenterMass[index] - mass) < Math.Abs(peakCenterMass[index + 1] - mass)) {
				return index;
			}
			return index + 1;
		}

		public Ms1CentroidList Extract(int[] indices) {
			Ms1CentroidList result = new Ms1CentroidList();
			result.peakCenterMass = ArrayUtil.SubArray(peakCenterMass, indices);
			result.peakMinMass = ArrayUtil.SubArray(peakMinMass, indices);
			result.peakMaxMass = ArrayUtil.SubArray(peakMaxMass, indices);
			result.peakIntensity = ArrayUtil.SubArray(peakIntensity, indices);
			result.peaks = ArrayUtil.SubArray(peaks, indices);
			return result;
		}

		internal void SetPeak(GrowablePeak peak, int peakIndex) {
			peaks[peakIndex] = peak;
		}

		public void Dispose() {
			peakCenterMass = null;
			peakMinMass = null;
			peakMaxMass = null;
			peakIntensity = null;
			for (int i = 0; i < peaks.Length; i++) {
				peaks[i] = null;
			}
			peaks = null;
		}
	}
}