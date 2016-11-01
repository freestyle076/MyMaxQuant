/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
namespace MaxQuant.Spec {
	public interface IPeakCollector {
		void AddPeak(Peak peak);
	}
}