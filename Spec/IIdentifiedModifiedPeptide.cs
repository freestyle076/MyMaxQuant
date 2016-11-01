/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using MaxQuant.Spec.Evidence;

namespace MaxQuant.Spec {
	public interface IIdentifiedModifiedPeptide {
		PeptideModificationCounts Counts { get; }
		SilacEvidence[] SilacEvidences { get; }
		IsotopeMsmsEvidence[] IsotopeMsmsEvidences { get; }
	}
}
