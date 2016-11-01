/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Collections.Generic;

namespace MaxQuant.Spec {
	public class GrowablePeak {
		private List<double> centerMz = new List<double>();
		private List<byte> massRange = new List<byte>();
		private List<float> maxMz = new List<float>();
		private List<float> minMz = new List<float>();
		private List<float> origIntensity = new List<float>();
		private List<int> scanIndex = new List<int>();

		public int LastScanIndex {
			get {
				if (centerMz.Count == 0) {
					return -1;
				}
				return scanIndex[centerMz.Count - 1];
			}
		}

		public void Add(int scanInd, int peakIndex, Ms1CentroidList peakList, byte massRang, double intensityNorm) {
			centerMz.Add(peakList.CenterMass(peakIndex));
			minMz.Add(peakList.MinMass(peakIndex));
			maxMz.Add(peakList.MaxMass(peakIndex));
			origIntensity.Add((float) (intensityNorm * peakList.GetIntensity(peakIndex)));
			scanIndex.Add(scanInd);
			massRange.Add(massRang);
			peakList.SetPeak(this, peakIndex);
		}

		public void Dispose() {
			centerMz = null;
			minMz = null;
			maxMz = null;
			origIntensity = null;
			scanIndex = null;
			massRange = null;
		}

		public bool IsDisposed() {
			return centerMz == null;
		}

		public Peak ToPeak() {
			return new Peak(centerMz.ToArray(), minMz.ToArray(), maxMz.ToArray(), new float[centerMz.Count],
			                origIntensity.ToArray(), scanIndex.ToArray(), massRange.ToArray());
		}
	}
}