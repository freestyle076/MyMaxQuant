/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.IO;
using MaxQuant.Spec;
using MaxQuant.Spec.Evidence;

namespace MaxQuant.Spec.Evidence {
	public class SilacMassEvidence : SilacEvidence {
		public SilacMassEvidence(BinaryReader reader) : base(reader) {
		}

		public SilacMassEvidence(int silacId, int fileIndex, SilacCluster silacCluster, int charge, double intensity0,
		                         double intensity1, double intensity2, double rt)
			: base(silacId, fileIndex, silacCluster, charge, intensity0,
			       intensity1, intensity2, rt) {
		}

		public override void Write(BinaryWriter writer) {
			writer.Write(silacPmfEvidence);
			base.Write(writer);
		}

		public override string Type {
			get { return "SILAC-MASS"; }
		}
	}
}