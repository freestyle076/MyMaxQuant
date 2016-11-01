/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.IO;

namespace MaxQuant.Spec {
	public class SilacCluster {
		private readonly int isotopeClusterIndex1 = -1;
		private readonly int isotopeClusterIndex2 = -1;
		private readonly int isotopeClusterIndex3 = -1;
		private readonly int offset10;
		private readonly int offset20;
		private readonly int offset21;
		private int[] aaCounts;
		private float elutionTime;
		private float elutionTimeWidth;
		private float intensity1 = float.NaN;
		private float intensity2 = float.NaN;
		private float intensity3 = float.NaN;
		private int isotopePatternStart1;
		private int isotopePatternStart2;
		private int isotopePatternStart3;
		private int[] labelCounts1;
		private int[] labelCounts2;
		private double mass = double.NaN;
		private double massDiff1;
		private double massDiff2;
		private double massError;
		private float normalizedRatio10 = float.NaN;
		private float normalizedRatio20 = float.NaN;
		private float normalizedRatio21 = float.NaN;
		private int npoints;
		private float ratio10 = float.NaN;
		private float ratio20 = float.NaN;
		private float ratio21 = float.NaN;
		private double ratioSignificanceA10;
		private double ratioSignificanceA20;
		private double ratioSignificanceA21;
		private double ratioSignificanceB10;
		private double ratioSignificanceB20;
		private double ratioSignificanceB21;
		private float silacIsotopeCorr;
		private float silacTimeCorr;
		private float sulphurCount;
		private float sulphurScore;
		private float theorIsotopeCorr;
		private double uncalibratedMass = double.NaN;

		public SilacCluster(SilacPair pair) {
			isotopeClusterIndex1 = pair.Index1;
			isotopeClusterIndex2 = pair.Index2;
			aaCounts = pair.LabelCounts;
			silacIsotopeCorr = (float) pair.IsotopeCorrelation;
			silacTimeCorr = (float) pair.TimeCorrelation;
			offset10 = pair.Offset;
		}

		public SilacCluster(SilacPair pair01, SilacPair pair02, SilacPair pair12) {
			if (pair01 != null) {
				isotopeClusterIndex1 = pair01.Index1;
				isotopeClusterIndex2 = pair01.Index2;
				aaCounts = pair01.LabelCounts;
				silacIsotopeCorr = (float) pair01.IsotopeCorrelation;
				silacTimeCorr = (float) pair01.TimeCorrelation;
				offset10 = pair01.Offset;
			}
			if (pair02 != null) {
				isotopeClusterIndex1 = pair02.Index1;
				isotopeClusterIndex3 = pair02.Index2;
				aaCounts = pair02.LabelCounts;
				silacIsotopeCorr = (float) pair02.IsotopeCorrelation;
				silacTimeCorr = (float) pair02.TimeCorrelation;
				offset20 = pair02.Offset;
			}
			if (pair12 != null) {
				isotopeClusterIndex2 = pair12.Index1;
				isotopeClusterIndex3 = pair12.Index2;
				aaCounts = pair12.LabelCounts;
				silacIsotopeCorr = (float) pair12.IsotopeCorrelation;
				silacTimeCorr = (float) pair12.TimeCorrelation;
				offset21 = pair12.Offset;
			}
		}

		public SilacCluster(int index1) {
			isotopeClusterIndex1 = index1;
			aaCounts = new int[0];
		}

		public SilacCluster(BinaryReader reader) {
			isotopeClusterIndex1 = reader.ReadInt32();
			isotopeClusterIndex2 = reader.ReadInt32();
			isotopeClusterIndex3 = reader.ReadInt32();
			int len = reader.ReadInt32();
			labelCounts1 = new int[len];
			for (int i = 0; i < labelCounts1.Length; i++) {
				labelCounts1[i] = reader.ReadInt32();
			}
			len = reader.ReadInt32();
			labelCounts2 = new int[len];
			for (int i = 0; i < labelCounts2.Length; i++) {
				labelCounts2[i] = reader.ReadInt32();
			}
			massDiff1 = reader.ReadDouble();
			massDiff2 = reader.ReadDouble();
			len = reader.ReadInt32();
			aaCounts = new int[len];
			for (int i = 0; i < aaCounts.Length; i++) {
				aaCounts[i] = reader.ReadInt32();
			}
			theorIsotopeCorr = reader.ReadSingle();
			silacIsotopeCorr = reader.ReadSingle();
			silacTimeCorr = reader.ReadSingle();
			isotopePatternStart1 = reader.ReadInt32();
			isotopePatternStart2 = reader.ReadInt32();
			isotopePatternStart3 = reader.ReadInt32();
			ratio10 = reader.ReadSingle();
			ratio20 = reader.ReadSingle();
			ratio21 = reader.ReadSingle();
			normalizedRatio10 = reader.ReadSingle();
			normalizedRatio20 = reader.ReadSingle();
			normalizedRatio21 = reader.ReadSingle();
			intensity1 = reader.ReadSingle();
			intensity2 = reader.ReadSingle();
			intensity3 = reader.ReadSingle();
			ratioSignificanceA10 = reader.ReadDouble();
			ratioSignificanceA20 = reader.ReadDouble();
			ratioSignificanceA21 = reader.ReadDouble();
			ratioSignificanceB10 = reader.ReadDouble();
			ratioSignificanceB20 = reader.ReadDouble();
			ratioSignificanceB21 = reader.ReadDouble();
			mass = reader.ReadDouble();
			uncalibratedMass = reader.ReadDouble();
			massError = reader.ReadDouble();
			sulphurScore = reader.ReadSingle();
			sulphurCount = reader.ReadSingle();
			elutionTime = reader.ReadSingle();
			offset10 = reader.ReadInt32();
			offset20 = reader.ReadInt32();
			offset21 = reader.ReadInt32();
			npoints = reader.ReadInt32();
			elutionTimeWidth = reader.ReadSingle();
		}

		public int[] AaCounts {
			get { return aaCounts; }
		}

		public double SulphurScore {
			set { sulphurScore = (float) value; }
			get { return sulphurScore; }
		}

		public double SulphurCount {
			set { sulphurCount = (float) value; }
			get { return sulphurCount; }
		}

		public int IsotopeClusterindex1 {
			get { return isotopeClusterIndex1; }
		}

		public int IsotopeClusterindex2 {
			get { return isotopeClusterIndex2; }
		}

		public int IsotopeClusterindex3 {
			get { return isotopeClusterIndex3; }
		}

		public int[] LabelCounts1 {
			get { return labelCounts1; }
		}

		public int[] LabelCounts2 {
			get { return labelCounts2; }
		}

		public int IsotopePatternStart1 {
			set { isotopePatternStart1 = value; }
			get { return isotopePatternStart1; }
		}

		public int IsotopePatternStart2 {
			set { isotopePatternStart2 = value; }
			get { return isotopePatternStart2; }
		}

		public int IsotopePatternStart3 {
			set { isotopePatternStart3 = value; }
			get { return isotopePatternStart3; }
		}

		public double Intensity1 {
			set { intensity1 = (float) value; }
			get { return intensity1; }
		}

		public double Intensity2 {
			set { intensity2 = (float) value; }
			get { return intensity2; }
		}

		public double Intensity3 {
			set { intensity3 = (float) value; }
			get { return intensity3; }
		}

		public double Ratio10 {
			set { ratio10 = (float) value; }
			get { return ratio10; }
		}

		public double Ratio20 {
			set { ratio20 = (float) value; }
			get { return ratio20; }
		}

		public double Ratio21 {
			set { ratio21 = (float) value; }
			get { return ratio21; }
		}

		public double NormalizedRatio10 {
			set { normalizedRatio10 = (float) value; }
			get { return normalizedRatio10; }
		}

		public double NormalizedRatio20 {
			set { normalizedRatio20 = (float) value; }
			get { return normalizedRatio20; }
		}

		public double NormalizedRatio21 {
			set { normalizedRatio21 = (float) value; }
			get { return normalizedRatio21; }
		}

		public double RatioSignificanceA10 {
			set { ratioSignificanceA10 = value; }
			get { return ratioSignificanceA10; }
		}

		public double RatioSignificanceA20 {
			set { ratioSignificanceA20 = value; }
			get { return ratioSignificanceA20; }
		}

		public double RatioSignificanceA21 {
			set { ratioSignificanceA21 = value; }
			get { return ratioSignificanceA21; }
		}

		public double RatioSignificanceB10 {
			set { ratioSignificanceB10 = value; }
			get { return ratioSignificanceB10; }
		}

		public double RatioSignificanceB20 {
			set { ratioSignificanceB20 = value; }
			get { return ratioSignificanceB20; }
		}

		public double RatioSignificanceB21 {
			set { ratioSignificanceB21 = value; }
			get { return ratioSignificanceB21; }
		}

		public double TheorIsotopeCorr {
			set { theorIsotopeCorr = (float) value; }
			get { return theorIsotopeCorr; }
		}

		public double SilacTimeCorr {
			set { silacTimeCorr = (float) value; }
			get { return silacTimeCorr; }
		}

		public double SilacIsotopeCorr {
			set { silacIsotopeCorr = (float) value; }
			get { return silacIsotopeCorr; }
		}

		public double MassDiff1 {
			set { massDiff1 = value; }
			get { return massDiff1; }
		}

		public double MassDiff2 {
			set { massDiff2 = value; }
			get { return massDiff2; }
		}

		public double ElutionTime {
			get { return elutionTime; }
			set { elutionTime = (float) value; }
		}

		public double ElutionTimeWidth {
			get { return elutionTimeWidth; }
			set { elutionTimeWidth = (float) value; }
		}

		public double Mass {
			get { return mass; }
			set { mass = value; }
		}

		public int Npoints {
			get { return npoints; }
			set { npoints = value; }
		}

		public double UncalibratedMass {
			get { return uncalibratedMass; }
			set { uncalibratedMass = value; }
		}

		public void Write(BinaryWriter writer) {
			writer.Write(isotopeClusterIndex1);
			writer.Write(isotopeClusterIndex2);
			writer.Write(isotopeClusterIndex3);
			int len = (labelCounts1 != null) ? labelCounts1.Length : 0;
			writer.Write(len);
			for (int i = 0; i < len; i++) {
				writer.Write(labelCounts1[i]);
			}
			len = (labelCounts2 != null) ? labelCounts2.Length : 0;
			writer.Write(len);
			for (int i = 0; i < len; i++) {
				writer.Write(labelCounts2[i]);
			}
			writer.Write(massDiff1);
			writer.Write(massDiff2);
			len = (aaCounts != null) ? aaCounts.Length : 0;
			writer.Write(len);
			for (int i = 0; i < len; i++) {
				writer.Write(aaCounts[i]);
			}
			writer.Write(theorIsotopeCorr);
			writer.Write(silacIsotopeCorr);
			writer.Write(silacTimeCorr);
			writer.Write(isotopePatternStart1);
			writer.Write(isotopePatternStart2);
			writer.Write(isotopePatternStart3);
			writer.Write(ratio10);
			writer.Write(ratio20);
			writer.Write(ratio21);
			writer.Write(normalizedRatio10);
			writer.Write(normalizedRatio20);
			writer.Write(normalizedRatio21);
			writer.Write(intensity1);
			writer.Write(intensity2);
			writer.Write(intensity3);
			writer.Write(ratioSignificanceA10);
			writer.Write(ratioSignificanceA20);
			writer.Write(ratioSignificanceA21);
			writer.Write(ratioSignificanceB10);
			writer.Write(ratioSignificanceB20);
			writer.Write(ratioSignificanceB21);
			writer.Write(mass);
			writer.Write(uncalibratedMass);
			writer.Write(massError);
			writer.Write(sulphurScore);
			writer.Write(sulphurCount);
			writer.Write(elutionTime);
			writer.Write(offset10);
			writer.Write(offset20);
			writer.Write(offset21);
			writer.Write(npoints);
			writer.Write(elutionTimeWidth);
		}

		public int[] GetIsotopeClusterIndices(SilacType type) {
			if (type == SilacType.Triplets) {
				return new int[] {isotopeClusterIndex1, isotopeClusterIndex2, isotopeClusterIndex3};
			}
			if (type == SilacType.Doublets) {
				return new int[] {isotopeClusterIndex1, isotopeClusterIndex2};
			}
			return new int[] {isotopeClusterIndex1};
		}

		public int[] GetIsotopePatternStarts(SilacType type) {
			if (type == SilacType.Triplets) {
				return new int[] {isotopePatternStart1, isotopePatternStart2, isotopePatternStart3};
			}
			if (type == SilacType.Doublets) {
				return new int[] {isotopePatternStart1, isotopePatternStart2};
			}
			return new int[] {isotopePatternStart1};
		}

		public double GetMassError(double factor) {
			return factor * massError;
		}

		public void SetMassError(double massError) {
			this.massError = massError;
		}

		public double[] GetMassDiffs(SilacType type) {
			if (type == SilacType.Triplets) {
				return new double[] {0, massDiff1, massDiff2};
			}
			if (type == SilacType.Doublets) {
				return new double[] {0, massDiff1};
			}
			return new double[] {0};
		}

		public int GetAaCount(int ind) {
			return aaCounts[ind];
		}

		public int[] GetAaCounts() {
			return aaCounts;
		}

		public int[] GetOffsets(SilacType type) {
			if (type == SilacType.Triplets) {
				return new int[] {offset10, offset20, offset21};
			}
			if (type == SilacType.Doublets) {
				return new int[] {offset10};
			}
			return new int[] {};
		}

		public void Dispose() {
			labelCounts1 = null;
			labelCounts2 = null;
			aaCounts = null;
		}

		public double[] GetIntensities(int silacCount) {
			if (silacCount == 1) {
				return new double[] {intensity1};
			}
			if (silacCount == 2) {
				return new double[] {intensity1, intensity2};
			}
			if (silacCount == 3) {
				return new double[] {intensity1, intensity2, intensity3};
			}
			throw new Exception("Invalid argument.");
		}
	}
}