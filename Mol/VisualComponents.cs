/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Drawing;

namespace MaxQuant.Mol
{
    public class VisualModification
    {
        private readonly int id = -1;
        private readonly int position;
        private readonly string label;
        private readonly string shortLabel;
        private readonly Type modificationTyp;

        public enum Type
        {
            known = 1,
            detected = 2
        } ;


        public VisualModification(string shortLabel, string label, int position, Type modificationTyp)
        {
            this.position = position;
            this.label = label;
            this.shortLabel = shortLabel;
            this.modificationTyp = modificationTyp;
        }

        public VisualModification(int id, string shortLabel, string label, int position, Type modificationTyp)
        {
            this.id = id;
            this.position = position;
            this.label = label;
            this.shortLabel = shortLabel;
            this.modificationTyp = modificationTyp;
        }

        public string Label
        {
            get { return label; }
        }

        public string ShortLabel
        {
            get { return shortLabel; }
        }

        public int Position
        {
            get { return position; }
        }

        public Type ModificationTyp
        {
            get { return modificationTyp; }
        }

        public int ID
        {
            get { return id; }
        }

        public override int GetHashCode()
        {
            int hash = id + position;
            hash += label.GetHashCode() + modificationTyp.GetHashCode();

            return hash;
        }

        public override bool Equals(object obj)
        {
            if (obj != null)
            {
                if (obj is VisualModification)
                {
                    VisualModification m = obj as VisualModification;
                    if (m.ID == this.ID && m.Position == this.Position)
                    {
                        if (m.Label.Equals(this.Label) && m.ModificationTyp.Equals(this.ModificationTyp))
                        {
                            return true;
                        }
                    }
                }
            }

            return false;
        }
    }

    public class VisualPeptide
    {
        private readonly Location locations;
        private string label;
        private readonly string shortLabel;
        private Color color;

        public VisualPeptide(string shortLabel, string label, Location locations)
        {
            this.locations = locations;
            this.label = label;
            this.shortLabel = shortLabel;
        }

        public string Label
        {
            get { return label; }
            set { label = value; }
        }

        public string ShortLabel
        {
            get { return shortLabel; }
        }

        public Location Locations
        {
            get { return locations; }
        }

        public Color Color
        {
            get { return color; }

            set { color = value; }
        }
    }

    public class VisualStructure
    {
        private readonly Location locations;
        private readonly string label;
        private readonly string shortLabel;
        private Color color;

        public VisualStructure(string shortLabel, string label, Location locations)
        {
            this.locations = locations;
            this.label = label;
            this.shortLabel = shortLabel;
        }

        public string Label
        {
            get { return label; }
        }

        public string ShortLabel
        {
            get { return shortLabel; }
        }

        public Location Locations
        {
            get { return locations; }
        }

        public Color Color
        {
            get { return color; }

            set { color = value; }
        }
    }

    public class VisualDomain
    {
        private readonly Location locations;
        private readonly string label;
        private readonly string shortLabel;
        private Color color;

        public VisualDomain(string shortLabel, string label, Location locations)
        {
            this.locations = locations;
            this.label = label;
            this.shortLabel = shortLabel;
        }

        public string Label
        {
            get { return label; }
        }

        public string ShortLabel
        {
            get { return shortLabel; }
        }

        public Location Locations
        {
            get { return locations; }
        }

        public Color Color
        {
            get { return color; }

            set { color = value; }
        }
    }
}