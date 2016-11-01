/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.IO;
using MaxQuant.Spec;
using MaxQuant.Spec.Evidence;

namespace MaxQuant.Spec.Evidence {
	public enum RatioType {
		Type10, Type20, Type21
	}

	public abstract class SilacEvidence : AbstractEvidence {
		private int silacId;
		private float ratio10;
		private float ratio20;
		private float ratio21;
		private float normalizedRatio10;
		private float normalizedRatio20;
		private float normalizedRatio21;
		private readonly float ratioSignificanceA10;
		private readonly float ratioSignificanceA20;
		private readonly float ratioSignificanceA21;
		private readonly float ratioSignificanceB10;
		private readonly float ratioSignificanceB20;
		private readonly float ratioSignificanceB21;
		private float intensity0;
		private float intensity1;
		private float intensity2;
		private readonly double calibratedMass;
		private readonly double uncalibratedMass;

		protected SilacEvidence(BinaryReader reader) : base(reader) {
			silacId = reader.ReadInt32();
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
			calibratedMass = reader.ReadDouble();
			uncalibratedMass = reader.ReadDouble();
		}

		protected SilacEvidence(int silacId, int rawFileIndex, SilacCluster silacCluster, int charge, 
		                        double intensity0, double intensity1, double intensity2, double calibRT)
			: base(rawFileIndex, (float)silacCluster.ElutionTime, (float)calibRT, charge) {
			this.silacId = silacId;
			ratio10 = (float)silacCluster.Ratio10;
			ratio20 = (float)silacCluster.Ratio20;
			ratio21 = (float)silacCluster.Ratio21;
			normalizedRatio10 = (float)silacCluster.NormalizedRatio10;
			normalizedRatio20 = (float)silacCluster.NormalizedRatio20;
			normalizedRatio21 = (float)silacCluster.NormalizedRatio21;
			ratioSignificanceA10 = (float)Math.Log(silacCluster.RatioSignificanceA10);
			ratioSignificanceA20 = (float)Math.Log(silacCluster.RatioSignificanceA20);
			ratioSignificanceA21 = (float)Math.Log(silacCluster.RatioSignificanceA21);
			ratioSignificanceB10 = (float)Math.Log(silacCluster.RatioSignificanceB10);
			ratioSignificanceB20 = (float)Math.Log(silacCluster.RatioSignificanceB20);
			ratioSignificanceB21 = (float)Math.Log(silacCluster.RatioSignificanceB21);
			this.intensity0 = (float)intensity0;
			this.intensity1 = (float)intensity1;
			this.intensity2 = (float)intensity2;
			calibratedMass = silacCluster.Mass;
			uncalibratedMass = silacCluster.UncalibratedMass;
		}

		public int SilacId {
			get { return silacId; }
			set { silacId = value; }
		}

		public double GetRatio(RatioType type) {
			switch(type) {
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
			set { ratio10 = (float)value; }
		}

		public double Ratio20 {
			get { return ratio20; }
			set { ratio20 = (float)value; }
		}

		public double Ratio21 {
			get { return ratio21; }
			set { ratio21 = (float)value; }
		}

		public double Intensity0 {
			get { return intensity0; }
			set { intensity0 = (float)value; }
		}

		public double Intensity1 {
			get { return intensity1; }
			set { intensity1 = (float)value; }
		}

		public double Intensity2 {
			get { return intensity2; }
			set { intensity2 = (float)value; }
		}

		public double GetIntensitySum(SilacType type, int index) {
			if(index >=0) {
				switch (index) {
					case 0:
						return intensity0;
					case 1:
						return intensity1;
					case 2:
						return intensity2;
					default:
						throw new Exception("Impossible.");
				}				
			}
			switch(type) {
				case SilacType.Singlets:
					return intensity0;
				case SilacType.Doublets:
					return intensity0 + intensity1;
				case SilacType.Triplets:
					return intensity0 + intensity1 + intensity2;
				default:
					throw new Exception("Impossible.");
			}
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
			set { normalizedRatio10 = (float)value; }
		}

		public double NormalizedRatio20 {
			get { return normalizedRatio20; }
			set { normalizedRatio20 = (float)value; }
		}

		public double NormalizedRatio21 {
			get { return normalizedRatio21; }
			set { normalizedRatio21 = (float)value; }
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

		public double CalibratedMass {
			get { return calibratedMass; }
		}

		public double UncalibratedMass {
			get { return uncalibratedMass; }
		}

		public override void Write(BinaryWriter writer) {
			base.Write(writer);
			writer.Write(silacId);
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
			writer.Write(calibratedMass);
			writer.Write(uncalibratedMass);
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
	}
}