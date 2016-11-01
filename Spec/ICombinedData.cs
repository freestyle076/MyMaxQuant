/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
namespace MaxQuant.Spec {
	public interface ICombinedData {
		int GetProteinGroupCount();
		IIdentifiedProteinGroup GetIProteinGroupAt(int i);
		int MinRatioCount { get; }
		RawFileInfo FileInfo { get; }
		QuantitationMode QuantitationMode { get; }
		bool UseCounterparts { get; }
		bool RestrictProteinQuantitation { get; }
		bool AreQuantifiableModifications(ushort[] types);
		IIdentifiedPeptide GetIPeptideAt(int i);
	}
}	
