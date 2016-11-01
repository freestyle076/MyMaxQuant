/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using MaxQuant.Mol;
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public interface IIdentifiedPeptide {
		void AddMascotPeptideHit(MascotPeptide[] p, int number, int index, MascotQueryType type, int id, int silacIndex,
		                         int isotopeId, SilacCluster cluster, IsotopeCluster isotopeCluster, double time,
		                         IPeakList list, double mz, double monotopicMz, ushort[] ushorts, double[] masses,
		                         float[] intensities, ReQuantitationResult result, bool quantify, HashSet<string> set,
		                         SilacType silacType, SilacLabel[] labels1, SilacLabel[] labels2, double tol, string unit,
		                         int topx, string[] mods);

		bool UniqueProtein {set; }
		bool UniqueGroup { get; set; }
		IIdentifiedModifiedPeptide[] IModifiedPeptides { get; }
	}
}