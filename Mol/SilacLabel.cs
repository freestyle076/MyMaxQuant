/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
namespace MaxQuant.Mol {
	public class SilacLabel {
		private readonly string mascotName;
		private readonly string shortName;

		public SilacLabel(string shortName, string mascotName) {
			this.shortName = shortName;
			this.mascotName = mascotName;
		}

		public string MascotName {
			get { return mascotName; }
		}

		public string ShortName {
			get { return shortName; }
		}

		public override bool Equals(object obj) {
			if (obj is SilacLabel) {
				SilacLabel other = (SilacLabel) obj;
				return other.shortName.Equals(shortName);
			}
			return false;
		}

		public override int GetHashCode() {
			return shortName.GetHashCode();
		}

		public override string ToString() {
			return shortName;
		}
	}
}