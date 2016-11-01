/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.IO;
using MaxQuant.Spec;

namespace MaxQuant.Mol {
	public class ModifiedPeptide {
		internal int peptideIndex;
		public PeptideModificationState modifications;
		public double mass;

		private ModifiedPeptide() {
		}

		public ModifiedPeptide(Peptide peptide, int peptideIndex) {
			this.peptideIndex = peptideIndex;
			modifications = new PeptideModificationState(peptide.Sequence.Length);
			mass = peptide.MonoIsotopicMass;
		}

		public ModifiedPeptide(DatabasePeptide peptide, int peptideIndex, string sequence) {
			this.peptideIndex = peptideIndex;
			modifications = new PeptideModificationState(sequence.Length);
			mass = peptide.GetMonoIsotopicMass(sequence);
		}

		public static ModifiedPeptide Read(BinaryReader reader) {
			try {
				ModifiedPeptide result = new ModifiedPeptide();
				result.peptideIndex = reader.ReadInt32();
				result.modifications = PeptideModificationState.Read(reader);
				result.mass = reader.ReadDouble();
				return result;
			} catch (EndOfStreamException) {
				return null;
			}
		}

		public void Write(BinaryWriter writer) {
			writer.Write(peptideIndex);
			modifications.Write(writer);
			writer.Write(mass);
		}

		internal ModifiedPeptide Clone() {
			ModifiedPeptide p = new ModifiedPeptide();
			p.peptideIndex = peptideIndex;
			p.modifications = modifications.Clone();
			p.mass = mass;
			return p;
		}
	}
}