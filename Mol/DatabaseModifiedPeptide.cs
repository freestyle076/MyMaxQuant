/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.IO;

namespace MaxQuant.Mol {
	public class DatabaseModifiedPeptide{
		public int peptideIndex;
		public double mass;

		private DatabaseModifiedPeptide() {
		}

		public DatabaseModifiedPeptide(int peptideIndex, double mass) {
			this.peptideIndex = peptideIndex;
			this.mass = mass;
		}

		public static DatabaseModifiedPeptide Read(BinaryReader reader) {
			try {
				DatabaseModifiedPeptide result = new DatabaseModifiedPeptide();
				result.peptideIndex = reader.ReadInt32();
				result.mass = reader.ReadDouble();
				return result;
			} catch (EndOfStreamException) {
				return null;
			}
		}

		public void Write(BinaryWriter writer) {
			writer.Write(peptideIndex);
			writer.Write(mass);
		}
	}
}