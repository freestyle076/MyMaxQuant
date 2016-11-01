/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using MaxQuant.Util;

namespace MaxQuant.Spec {
	public interface IIdentifiedProteinGroup {
		void GetRatio10(ICombinedData data, out double ratio, out double normalizedRatio, out int count, out double dev, out int outliers, out double pval);
		double GetIntensity0(ICombinedData data);
		double GetIntensity1(ICombinedData data);
		double GetIntensity2(ICombinedData data);
		double RatioSignificanceA10 { get; set; }
		double RatioSignificanceB10 { get; set; }
		double RatioSignificanceA20 { get; set; }
		double RatioSignificanceB20 { get; set; }
		double RatioSignificanceA21 { get; set; }
		double RatioSignificanceB21 { get; set; }
		void AddRatioSignificanceA10(string exp, double d);
		void AddRatioSignificanceB10(string exp, double d);
		void GetRatio10(ICombinedData data, out double ratio, out double normalizedRatio, out int count, out double dev, out int outliers, HashSet<int> set, out double pval);
		void GetRatio20(ICombinedData data, out double ratio, out double normalizedRatio, out int count, out double dev, out int outliers, out double pval);
		void GetRatio20(ICombinedData data, out double ratio, out double normalizedRatio, out int count, out double dev, out int outliers, HashSet<int> set, out double pval);
		void AddRatioSignificanceA20(string exp, double d);
		void AddRatioSignificanceB20(string exp, double d);
		void GetRatio21(ICombinedData data, out double ratio, out double normalizedRatio, out int count, out double dev, out int outliers, out double pval);
		void AddRatioSignificanceA21(string exp, double d);
		void AddRatioSignificanceB21(string exp, double d);
		void GetRatio21(ICombinedData data, out double ratio, out double normalizedRatio, out int count, out double dev, out int outliers, HashSet<int> set, out double pval);
	}
}
