/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Text;

namespace MaxQuant.Mol
{
    class SecondaryStructureElement
    {
        public static readonly SecondaryStructureElement Helix = new SecondaryStructureElement("helix", 'H');
        public static readonly SecondaryStructureElement Sheet = new SecondaryStructureElement("sheet", 'E');
        public static readonly SecondaryStructureElement RandomCoil = new SecondaryStructureElement("random coil", 'C');

        public static readonly SecondaryStructureElement[] sstructureElements = new SecondaryStructureElement[] { Helix, Sheet, RandomCoil };

        private readonly string name;
        private readonly char letter;

        public static readonly string singleLetter = ExtractSingleLetter();

        internal SecondaryStructureElement(string name, char letter)
        {
            this.name = name;
            this.letter = letter;
        }

        public char Letter
        {
            get { return letter; }
        }

        public string Name
        {
            get { return name; }
        }

        private static string ExtractSingleLetter()
        {
            StringBuilder result = new StringBuilder();
            for (int i = 0; i < sstructureElements.Length; i++)
            {
                result.Append(sstructureElements[i].letter);
            }
            return result.ToString();
        }
    }
}
