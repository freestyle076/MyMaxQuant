/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.IO;
using MaxQuant.Spec;
using MaxQuant.Spec.Evidence;

namespace MaxQuant.Spec.Evidence {
	public class SilacMsmsEvidence : SilacEvidence {
		public SilacMsmsEvidence(BinaryReader reader) : base(reader) {
		}

		public SilacMsmsEvidence(int silacId, int rawFileIndex, SilacCluster silacCluster, int charge, double intensity0,
		                         double intensity1, double intensity2, double calibRT)
			:
				base(silacId, rawFileIndex, silacCluster, charge, intensity0, intensity1, intensity2, calibRT) {
		}

		public override void Write(BinaryWriter writer) {
			writer.Write(silacMsmsEvidence);
			base.Write(writer);
		}

		public override string Type {
			get { return "SILAC-MSMS"; }
		}
	}
}