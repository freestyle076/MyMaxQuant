/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.IO;
using MaxQuant.Spec;
using MaxQuant.Spec.Evidence;

namespace MaxQuant.Spec.Evidence {
	public class IsotopeMsmsEvidence : AbstractEvidence {
		private readonly int isotopeId;
		private float ratio10 = float.NaN;
		private float ratio20 = float.NaN;
		private float ratio21 = float.NaN;
		private float normalizedRatio10 = float.NaN;
		private float normalizedRatio20 = float.NaN;
		private float normalizedRatio21 = float.NaN;
		private readonly float ratioSignificanceA10 = float.NaN;
		private readonly float ratioSignificanceA20 = float.NaN;
		private readonly float ratioSignificanceA21 = float.NaN;
		private readonly float ratioSignificanceB10 = float.NaN;
		private readonly float ratioSignificanceB20 = float.NaN;
		private readonly float ratioSignificanceB21 = float.NaN;
		private float intensity0 = float.NaN;
		private float intensity1 = float.NaN;
		private float intensity2 = float.NaN;
		private readonly float intensity;

		public IsotopeMsmsEvidence(BinaryReader reader) : base(reader) {
			isotopeId = reader.ReadInt32();
			ratio10 = reader.ReadSingle();
			ratio20 = reader.ReadSingle();
			ratio21 = reader.ReadSingle();
			normalizedRatio10 = reader.ReadSingle();
			normalizedRatio20 = reader.ReadSingle();
			normalizedRatio21 = reader.ReadSingle();
			ratioSignificanceA10 = reader.ReadSingle();
			ratioSignificanceA20 = reader.ReadSingle();
			ratioSignificanceA21 = reader.ReadSingle();
			ratioSignificanceB10 = reader.ReadSingle();
			ratioSignificanceB20 = reader.ReadSingle();
			ratioSignificanceB21 = reader.ReadSingle();
			intensity0 = reader.ReadSingle();
			intensity1 = reader.ReadSingle();
			intensity2 = reader.ReadSingle();
			intensity = reader.ReadSingle();
		}

		public IsotopeMsmsEvidence(int isotopeId, int fileIndex, IsotopeCluster isotopeCluster,
		                           double intensity, double calibRT, ReQuantitationResult reQuantitationResult,
		                           SilacType type, bool reQuantify)
			: base(
				fileIndex, (isotopeCluster != null) ? (float) isotopeCluster.ElutionTime : float.NaN, (float) calibRT,
				(isotopeCluster != null) ? isotopeCluster.Charge : 0) {
			this.isotopeId = isotopeId;
			this.intensity = (float) intensity;
			if (reQuantify && reQuantitationResult.ContainsIsotopeIndex(isotopeId, type)) {
				if (type == SilacType.Doublets) {
					ratio10 = (float) reQuantitationResult.GetDoubletRatio(isotopeId);
					normalizedRatio10 = (float) reQuantitationResult.GetDoubletNormalizedRatio(isotopeId);
					ratioSignificanceA10 = (float) Math.Log(reQuantitationResult.GetDoubletRatioSignificanceA(isotopeId));
					ratioSignificanceB10 = (float) Math.Log(reQuantitationResult.GetDoubletRatioSignificanceB(isotopeId));
					double[] ii = reQuantitationResult.GetDoubletIntensities(isotopeId);
					intensity0 = (float) ii[0];
					intensity1 = (float) ii[1];
				}
				if (type == SilacType.Triplets) {
					double[] r = reQuantitationResult.GetTripletRatio(isotopeId);
					ratio10 = (float) r[0];
					ratio20 = (float) r[1];
					ratio21 = (float) r[2];
					double[] rn = reQuantitationResult.GetTripletNormalizedRatio(isotopeId);
					normalizedRatio10 = (float) rn[0];
					normalizedRatio20 = (float) rn[1];
					normalizedRatio21 = (float) rn[2];
					double[] rsa = reQuantitationResult.GetTripletRatioSignificanceA(isotopeId);
					ratioSignificanceA10 = (float) Math.Log(rsa[0]);
					ratioSignificanceA20 = (float) Math.Log(rsa[1]);
					ratioSignificanceA21 = (float) Math.Log(rsa[2]);
					double[] rsb = reQuantitationResult.GetTripletRatioSignificanceB(isotopeId);
					ratioSignificanceB10 = (float) Math.Log(rsb[0]);
					ratioSignificanceB20 = (float) Math.Log(rsb[1]);
					ratioSignificanceB21 = (float) Math.Log(rsb[2]);
					double[] ii = reQuantitationResult.GetTripletIntensities(isotopeId);
					intensity0 = (float) ii[0];
					intensity1 = (float) ii[1];
					intensity2 = (float) ii[2];
				}
			}
		}

		public int IsotopeId {
			get { return isotopeId; }
		}

		public double Intensity {
			get { return intensity; }
		}

		public override string Type {
			get { return "ISO-MSMS"; }
		}

		public override void Write(BinaryWriter writer) {
			writer.Write(isotopeMsmsEvidence);
			base.Write(writer);
			writer.Write(isotopeId);
			writer.Write(ratio10);
			writer.Write(ratio20);
			writer.Write(ratio21);
			writer.Write(normalizedRatio10);
			writer.Write(normalizedRatio20);
			writer.Write(normalizedRatio21);
			writer.Write(ratioSignificanceA10);
			writer.Write(ratioSignificanceA20);
			writer.Write(ratioSignificanceA21);
			writer.Write(ratioSignificanceB10);
			writer.Write(ratioSignificanceB20);
			writer.Write(ratioSignificanceB21);
			writer.Write(intensity0);
			writer.Write(intensity1);
			writer.Write(intensity2);
			writer.Write(intensity);
		}

		public double GetRatio(RatioType type) {
			switch (type) {
				case RatioType.Type10:
					return ratio10;
				case RatioType.Type20:
					return ratio20;
				case RatioType.Type21:
					return ratio21;
			}
			throw new Exception("Never get here.");
		}

		public double Ratio10 {
			get { return ratio10; }
			set { ratio10 = (float) value; }
		}

		public double Ratio20 {
			get { return ratio20; }
			set { ratio20 = (float) value; }
		}

		public double Ratio21 {
			get { return ratio21; }
			set { ratio21 = (float) value; }
		}

		public double Intensity0 {
			get { return intensity0; }
			set { intensity0 = (float) value; }
		}

		public double Intensity1 {
			get { return intensity1; }
			set { intensity1 = (float) value; }
		}

		public double Intensity2 {
			get { return intensity2; }
			set { intensity2 = (float) value; }
		}

		public double GetNormalizedRatio(RatioType type) {
			switch (type) {
				case RatioType.Type10:
					return normalizedRatio10;
				case RatioType.Type20:
					return normalizedRatio20;
				case RatioType.Type21:
					return normalizedRatio21;
			}
			throw new Exception("Never get here.");
		}

		public double NormalizedRatio10 {
			get { return normalizedRatio10; }
			set { normalizedRatio10 = (float) value; }
		}

		public double NormalizedRatio20 {
			get { return normalizedRatio20; }
			set { normalizedRatio20 = (float) value; }
		}

		public double NormalizedRatio21 {
			get { return normalizedRatio21; }
			set { normalizedRatio21 = (float) value; }
		}

		public double RatioSignificanceA10 {
			get { return Math.Exp(ratioSignificanceA10); }
		}

		public double RatioSignificanceA20 {
			get { return Math.Exp(ratioSignificanceA20); }
		}

		public double RatioSignificanceA21 {
			get { return Math.Exp(ratioSignificanceA21); }
		}

		public double RatioSignificanceB10 {
			get { return Math.Exp(ratioSignificanceB10); }
		}

		public double RatioSignificanceB20 {
			get { return Math.Exp(ratioSignificanceB20); }
		}

		public double RatioSignificanceB21 {
			get { return Math.Exp(ratioSignificanceB21); }
		}

		public double GetIntensity(int which) {
			switch (which) {
				case 0:
					return intensity0;
				case 1:
					return intensity1;
				case 2:
					return intensity2;
			}
			throw new Exception("Never get here.");
		}

		public bool IsRequantified() {
			return !double.IsNaN(ratio10);
		}
	}
}