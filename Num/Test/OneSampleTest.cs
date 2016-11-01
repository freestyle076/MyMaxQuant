/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
namespace MaxQuant.Num.Test {
	public abstract class OneSampleTest {
		public abstract void Test(double[] data, double mean, out double statistic, out double bothTails,
		                          out double leftTail, out double rightTail);
	}
}