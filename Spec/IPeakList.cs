/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Collections.Generic;

namespace MaxQuant.Spec {
	public interface IPeakList {
		Peak GetPeakKeep(int i);
		Peak GetPeakDiscard(int i);
		double GetMzDiscard(int ind, double[,] par, double[,] calibrationPar, out int np);
		double GetMz(int ind, double[,] par, double[,] calibrationPar, out int np);
		int GetSilacCharge(int id);
		double GetSilacMass(int id, int index, SilacType type);
		int[] GetSilacInfoForMsmsScanNumber(int number);
		int GetIsotopeIndexForMsmsScanNumber(int number);
		IsotopeCluster GetIsotopeCluster(int index);
		double GetMassEstimate(int[] ints, int[] ints1, double[] doubles, bool b);
		int GetPeakIndexForMsmsScanNumber(int number);
		double GetMs2Mz(int i);
		double GetMz(int index);
		double GetIsotopeClusterIntensity(int ind);

		Dictionary<int, double>[] Vectorize(IsotopeCluster c, double[][] d1, double[] masses1, float[][] intens1,
		                                    int[][] indices1, double[][] d2);

		Dictionary<int, double>[] Vectorize(double[] masses1, float[][] intens1, int[][] indices1, double[][] d2,
		                                    double[] masses2, float[][] intens2, int[][] indices2, double[][] d3, int charge);

		double[] MS2MonoisotopicMz { get; }
		int GetMs2IndexFromScanNumber(int number);
		double GetMs2Rt(int ms2ind);
		SilacCluster GetSilacCluster(int id);
		double GetCalibratedRT(double time);
	}
}