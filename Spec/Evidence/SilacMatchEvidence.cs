/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.IO;
using MaxQuant.Spec;
using MaxQuant.Spec.Evidence;

namespace MaxQuant.Spec.Evidence {
	public class SilacMatchEvidence : SilacEvidence {

		private readonly double massErrPpm;
		
		public SilacMatchEvidence(BinaryReader reader) : base(reader) {
			massErrPpm = reader.ReadDouble();
		}

		public SilacMatchEvidence(int silacId, int rawFileIndex, SilacCluster silacCluster,
		                          int charge, double intensity0, double intensity1,
		                          double intensity2, double calibRT, MascotPeptide peptide)
			:
				base(silacId, rawFileIndex, silacCluster, charge, intensity0,
				     intensity1, intensity2, calibRT) {
			massErrPpm = (silacCluster.Mass - peptide.Mass) / peptide.Mass * 1e6;
		}

		public override void Write(BinaryWriter writer) {
			writer.Write(SilacMatchEvidence);
			base.Write(writer);
			writer.Write(massErrPpm);
		}

		public override string Type {
			get { return "SILAC-MATCH"; }
		}

		public override double GetMassErrorPpm() {
			return massErrPpm;
		}
	}
}