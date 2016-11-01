/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.IO;
using MaxQuant.Spec.Evidence;

namespace MaxQuant.Spec.Evidence {
	public class PeakMsmsEvidence : AbstractEvidence {

		public PeakMsmsEvidence(BinaryReader reader) : base(reader) {
		}

		public PeakMsmsEvidence(int fileIndex, int charge, double elutionTime, double calibratedElutionTime)
			: base(fileIndex, (float)elutionTime, (float)calibratedElutionTime, charge) {
		}
		public override void Write(BinaryWriter writer) {
			writer.Write(peakMsmsEvidence);
			base.Write(writer);
		}
		public override string Type {
			get { return "MSMS"; }
		}
	}
}