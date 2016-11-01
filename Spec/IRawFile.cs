/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
namespace MaxQuant.Spec {
	public interface IRawFile {
		int MS1Count { get; }
		int NumberOfMS1MassRanges { get; }
		int MS2Count { get; }
		double[] GetMS1MassRange(int i);
		byte GetMS1MassRangeIndex(int i);
		double[] GetMS1TimeSpan(int i);
		Spectrum GetMS1Spectrum(int j, bool subtractBackground, int quantile);
		Spectrum GetMS2Spectrum(int i);
		int GetScanNumberFromMs2Index(int i);
		double GetMS2MonoisotopicMz(int i);
		SignalType GetMS2SignalType(int i);
	}
}