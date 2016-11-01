/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
namespace MaxQuant.Spec {
	public class SilacPair {
		private readonly int isoCluster1;
		private readonly int isoCluster2;

		private readonly double isotopeCorrelation;
		private readonly int[] labelCounts;
		private readonly int offset;
		private readonly double ratio;
		private readonly double timeCorrelation;

		public SilacPair(int isoCluster1, int isoCluster2, int[] labelCounts, double isotopeCorrelation,
		                 double timeCorrelation, double ratio, int offset) {
			this.isoCluster1 = isoCluster1;
			this.isoCluster2 = isoCluster2;
			this.labelCounts = labelCounts;
			this.isotopeCorrelation = isotopeCorrelation;
			this.timeCorrelation = timeCorrelation;
			this.ratio = ratio;
			this.offset = offset;
		}

		public int Index1 {
			get { return isoCluster1; }
		}

		public int Index2 {
			get { return isoCluster2; }
		}

		public int[] LabelCounts {
			get { return labelCounts; }
		}

		public double TimeCorrelation {
			get { return timeCorrelation; }
		}

		public double Ratio {
			get { return ratio; }
		}

		public double IsotopeCorrelation {
			get { return isotopeCorrelation; }
		}

		public int Offset {
			get { return offset; }
		}

		public bool EqualLabeledAaCounts(SilacPair other) {
			for (int i = 0; i < labelCounts.Length; i++) {
				if (labelCounts[i] != other.LabelCounts[i]) {
					return false;
				}
			}
			return true;
		}
	}
}