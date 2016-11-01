/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.IO;

namespace MaxQuant.Spec {
	public class IsotopeCluster {
		private readonly int charge;
		private int contaminantIndex = -1;

		private float elutionTime;
		private float elutionTimeWidth;
		private float intensity;
		private float isotopeCorrelation;
		private int isotopePatternStart;
		private double massError;
		private int[] members;
		private int polymerIndex = -1;
		private int silacClusterIndex = -1;

		public IsotopeCluster(BinaryReader reader) {
			int count = reader.ReadInt32();
			members = new int[count];
			for (int i = 0; i < count; i++) {
				members[i] = reader.ReadInt32();
			}
			charge = reader.ReadInt32();
			silacClusterIndex = reader.ReadInt32();
			polymerIndex = reader.ReadInt32();
			contaminantIndex = reader.ReadInt32();
			isotopePatternStart = reader.ReadInt32();
			massError = reader.ReadDouble();
			isotopeCorrelation = reader.ReadSingle();
			elutionTime = reader.ReadSingle();
			intensity = reader.ReadSingle();
			elutionTimeWidth = reader.ReadSingle();
		}

		public IsotopeCluster(int[] members, int charge) {
			this.members = members;
			this.charge = charge;
		}

		public int Count {
			get { return members.Length; }
		}

		public int[] Members {
			get { return members; }
		}

		public int Charge {
			get { return charge; }
		}

		public int SilacCluster {
			get { return silacClusterIndex; }
			set { silacClusterIndex = value; }
		}

		public int IsotopePatternStart {
			get { return isotopePatternStart; }
			set { isotopePatternStart = value; }
		}

		public double IsotopeCorrelation {
			get { return isotopeCorrelation; }
			set { isotopeCorrelation = (float) value; }
		}

		public int PolymerIndex {
			get { return polymerIndex; }
			set { polymerIndex = value; }
		}

		public int ContaminantIndex {
			get { return contaminantIndex; }
			set { contaminantIndex = value; }
		}

		public double ElutionTime {
			get { return elutionTime; }
			set { elutionTime = (float) value; }
		}

		public double ElutionTimeWidth {
			get { return elutionTimeWidth; }
			set { elutionTimeWidth = (float) value; }
		}

		public double Intensity {
			get { return intensity; }
			set { intensity = (float) value; }
		}

		public void Write(BinaryWriter writer) {
			writer.Write(members.Length);
			for (int i = 0; i < members.Length; i++) {
				writer.Write(members[i]);
			}
			writer.Write(charge);
			writer.Write(silacClusterIndex);
			writer.Write(polymerIndex);
			writer.Write(contaminantIndex);
			writer.Write(isotopePatternStart);
			writer.Write(massError);
			writer.Write(isotopeCorrelation);
			writer.Write(elutionTime);
			writer.Write(intensity);
			writer.Write(elutionTimeWidth);
		}

		public double GetMassError(double factor) {
			return factor * massError;
		}

		public void SetMassError(double massErr) {
			massError = massErr;
		}

		public void Dispose() {
			members = null;
		}
	}
}