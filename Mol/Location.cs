/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
namespace MaxQuant.Mol
{
    public class Location
    {
        private readonly int start;
        private readonly int length;
        private readonly int end;

        public Location(int start, int length)
        {
            this.start = start;
            this.length = length;
            end = start + length - 1;            
        }

        public int Start
        {
            get
            {
                return start;
            }
        }

        public int End
        {
            get
            {
                return end;
            }
        }

        public int Length
        {
            get
            {
                return length;
            }
        }
    }
}
