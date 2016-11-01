/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Collections.Generic;

namespace MaxQuant.Mol
{
    public class Sequence
    {
        private readonly string proteinSequence;
        private List<VisualModification> modifications;
        private List<VisualStructure> structures;
        private List<VisualDomain> domains;
        private List<VisualPeptide> identifiedPeptides;

        public Sequence(string sequence)
        {
            proteinSequence = sequence;
        }

        public string ProteinSequence
        {
            get
            {
                return proteinSequence;
            }
        }

        public List<VisualModification> Modifications
        {
            get
            {
                if (modifications == null)
                {
                    modifications = new List<VisualModification>();
                }
                return modifications;
            }

            set
            {
                modifications = value;
            }
        }

        public List<VisualPeptide> IdentifiedPeptides
        {
            get
            {
                if (identifiedPeptides == null)
                {
                    identifiedPeptides = new List<VisualPeptide>();
                }
                return identifiedPeptides;
            }

            set
            {
                identifiedPeptides = value;
            }
        }

        public List<VisualStructure> Structures
        {
            get
            {
                if (structures== null)
                {
                    structures = new List<VisualStructure>();
                }
                return structures;
            }

            set
            {
                structures = value;
            }
        }

        public List<VisualDomain> Domain
        {
            get
            {
                if (domains == null)
                {
                    domains = new List<VisualDomain>();
                }
                return domains;
            }

            set
            {
                domains = value;
            }
        }

        public int Length
        {
            get
            {
                return proteinSequence.Length;
            }
        }
    }
}
