/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.IO;
using MaxQuant.Util;

namespace MaxQuant.Mol
{
    public class Protein
    {
        private string accession;
        private string description;
        private string sequence;
        private Sequence sequenceData;
        private readonly double molecularWeight;

        public Protein(string sequence, string accession, string description)
        {
            this.sequence = sequence;
            this.accession = accession;
            this.description = description;
            molecularWeight = AminoAcid.CalcMolecularWeight(sequence);
        }

        public Protein(BinaryReader reader)
        {
            sequence = AASequence.Read(reader).ToString();
            sequenceData = new Sequence(sequence);
            accession = FileUtil.ReadString(reader);
            description = FileUtil.ReadString(reader);
            molecularWeight = reader.ReadDouble();
        }

        public void Write(BinaryWriter writer)
        {
            new AASequence(sequence).Write(writer);
            FileUtil.WriteString(accession, writer);
            FileUtil.WriteString(description, writer);
            writer.Write(molecularWeight);
        }

        public string Sequence
        {
            get { return sequence; }
            set { sequence = value; }
        }

        public string Accession
        {
            get { return accession; }
            set { accession = value; }
        }

        public string Description
        {
            get { return description; }
            set { description = value; }
        }

        public double MolecularWeight
        {
            get { return molecularWeight; }
        }

        public Sequence SequenceData
        {
            get { return sequenceData; }
            set { sequenceData = value; }
        }
    }
}
