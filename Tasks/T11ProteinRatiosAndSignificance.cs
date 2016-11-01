/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using MaxQuant.Num;
using MaxQuant.Num.Test;
using MaxQuant.Spec;
using MaxQuant.Spec.Evidence;
using MaxQuant.Util;

namespace MaxQuant.Tasks {
	public static class T11ProteinRatiosAndSignificance {
		public static void CalcRatioSignificanceProteinGroups(SilacType type, ICombinedData combinedData, RawFileInfo fileInfo) {
			if (type == SilacType.Doublets) {
				{
					List<double> lNormRatio = new List<double>();
					List<double> lRatio = new List<double>();
					List<double> lIntensity = new List<double>();
					List<int> indices = new List<int>();
					for (int i = 0; i < combinedData.GetProteinGroupCount(); i++) {
						IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(i);
						double ratio;
						double normalizedRatio;
						int count;
						double dev;
						int nOutliers;
						double pval;
						group.GetRatio10(combinedData, out ratio, out normalizedRatio, out count, out dev, out nOutliers, out pval);
						double intensity = group.GetIntensity0(combinedData) + group.GetIntensity1(combinedData);
						if (normalizedRatio > 0 && !double.IsInfinity(normalizedRatio) && !double.IsNaN(normalizedRatio) && intensity > 0) {
							lNormRatio.Add(Math.Log(normalizedRatio));
							lRatio.Add(Math.Log(ratio));
							lIntensity.Add(Math.Log(intensity));
							indices.Add(i);
						}
					}
					double[] ratioSignificanceA = NumUtil.MovingBoxPlot(lNormRatio.ToArray(), lIntensity.ToArray(), 1);
					double[] ratioSignificanceB = NumUtil.MovingBoxPlot(lNormRatio.ToArray(), lIntensity.ToArray(), -1);
					for (int i = 0; i < indices.Count; i++) {
						IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(indices[i]);
						group.RatioSignificanceA10 = ratioSignificanceA[i];
						group.RatioSignificanceB10 = ratioSignificanceB[i];
					}
				}
				if (fileInfo.HasExperiment) {
					string[] exps = fileInfo.GetAllExperimentValues();
					foreach (string exp in exps) {
						List<double> lratio = new List<double>();
						List<double> lIntensity = new List<double>();
						List<int> indices = new List<int>();
						for (int i = 0; i < combinedData.GetProteinGroupCount(); i++) {
							IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(i);
							double ratio;
							double normalizedRatio;
							int count;
							double dev;
							int nOutliers;
							double pval;
							group.GetRatio10(combinedData, out ratio, out normalizedRatio, out count, out dev, out nOutliers,
											 fileInfo.GetRawFileIndicesFromExperiment(exp), out pval);
							double intensity = group.GetIntensity0(combinedData) + group.GetIntensity1(combinedData);
							if (normalizedRatio > 0 && !double.IsInfinity(normalizedRatio) && !double.IsNaN(normalizedRatio) && intensity > 0) {
								lratio.Add(Math.Log(normalizedRatio));
								lIntensity.Add(Math.Log(intensity));
								indices.Add(i);
							}
						}
						double[] ratioSignificanceA = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), 1);
						double[] ratioSignificanceB = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), -1);
						for (int i = 0; i < indices.Count; i++) {
							IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(indices[i]);
							group.AddRatioSignificanceA10(exp, ratioSignificanceA[i]);
							group.AddRatioSignificanceB10(exp, ratioSignificanceB[i]);
						}
					}
				}
			}
			if (type == SilacType.Triplets) {
				{
					List<double> lratio = new List<double>();
					List<double> lIntensity = new List<double>();
					List<int> indices = new List<int>();
					for (int i = 0; i < combinedData.GetProteinGroupCount(); i++) {
						IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(i);
						double ratio;
						double normalizedRatio;
						int count;
						double dev;
						int nOutliers;
						double pval;
						group.GetRatio10(combinedData, out ratio, out normalizedRatio, out count, out dev, out nOutliers, out pval);
						double intensity = group.GetIntensity0(combinedData) + group.GetIntensity1(combinedData) +
										   group.GetIntensity2(combinedData);
						if (normalizedRatio > 0 && !double.IsInfinity(normalizedRatio) && !double.IsNaN(normalizedRatio) && intensity > 0) {
							lratio.Add(Math.Log(normalizedRatio));
							lIntensity.Add(Math.Log(intensity));
							indices.Add(i);
						}
					}
					double[] ratioSignificanceA = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), 1);
					double[] ratioSignificanceB = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), -1);
					for (int i = 0; i < indices.Count; i++) {
						IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(indices[i]);
						group.RatioSignificanceA10 = ratioSignificanceA[i];
						group.RatioSignificanceB10 = ratioSignificanceB[i];
					}
				}
				if (fileInfo.HasExperiment) {
					string[] exps = fileInfo.GetAllExperimentValues();
					foreach (string exp in exps) {
						List<double> lratio = new List<double>();
						List<double> lIntensity = new List<double>();
						List<int> indices = new List<int>();
						for (int i = 0; i < combinedData.GetProteinGroupCount(); i++) {
							IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(i);
							double ratio;
							double normalizedRatio;
							int count;
							double dev;
							int nOutliers;
							double pval;
							group.GetRatio10(combinedData, out ratio, out normalizedRatio, out count, out dev, out nOutliers,
											 fileInfo.GetRawFileIndicesFromExperiment(exp), out pval);
							double intensity = group.GetIntensity0(combinedData) + group.GetIntensity1(combinedData) +
											   group.GetIntensity2(combinedData);
							if (normalizedRatio > 0 && !double.IsInfinity(normalizedRatio) && !double.IsNaN(normalizedRatio) && intensity > 0) {
								lratio.Add(Math.Log(normalizedRatio));
								lIntensity.Add(Math.Log(intensity));
								indices.Add(i);
							}
						}
						double[] ratioSignificanceA = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), 1);
						double[] ratioSignificanceB = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), -1);
						for (int i = 0; i < indices.Count; i++) {
							IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(indices[i]);
							group.AddRatioSignificanceA10(exp, ratioSignificanceA[i]);
							group.AddRatioSignificanceB10(exp, ratioSignificanceB[i]);
						}
					}
				}
				{
					List<double> lratio = new List<double>();
					List<double> lIntensity = new List<double>();
					List<int> indices = new List<int>();
					for (int i = 0; i < combinedData.GetProteinGroupCount(); i++) {
						IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(i);
						double ratio;
						double normalizedRatio;
						int count;
						double dev;
						int nOutliers;
						double pval;
						group.GetRatio20(combinedData, out ratio, out normalizedRatio, out count, out dev, out nOutliers, out pval);
						double intensity = group.GetIntensity0(combinedData) + group.GetIntensity1(combinedData) +
										   group.GetIntensity2(combinedData);
						if (normalizedRatio > 0 && !double.IsInfinity(normalizedRatio) && !double.IsNaN(normalizedRatio) && intensity > 0) {
							lratio.Add(Math.Log(normalizedRatio));
							lIntensity.Add(Math.Log(intensity));
							indices.Add(i);
						}
					}
					double[] ratioSignificanceA = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), 1);
					double[] ratioSignificanceB = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), -1);
					for (int i = 0; i < indices.Count; i++) {
						IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(indices[i]);
						group.RatioSignificanceA20 = ratioSignificanceA[i];
						group.RatioSignificanceB20 = ratioSignificanceB[i];
					}
				}
				if (fileInfo.HasExperiment) {
					string[] exps = fileInfo.GetAllExperimentValues();
					foreach (string exp in exps) {
						List<double> lratio = new List<double>();
						List<double> lIntensity = new List<double>();
						List<int> indices = new List<int>();
						for (int i = 0; i < combinedData.GetProteinGroupCount(); i++) {
							IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(i);
							double ratio;
							double normalizedRatio;
							int count;
							double dev;
							int nOutliers;
							double pval;
							group.GetRatio20(combinedData, out ratio, out normalizedRatio, out count, out dev, out nOutliers,
											 fileInfo.GetRawFileIndicesFromExperiment(exp), out pval);
							double intensity = group.GetIntensity0(combinedData) + group.GetIntensity1(combinedData) +
											   group.GetIntensity2(combinedData);
							if (normalizedRatio > 0 && !double.IsInfinity(normalizedRatio) && !double.IsNaN(normalizedRatio) && intensity > 0) {
								lratio.Add(Math.Log(normalizedRatio));
								lIntensity.Add(Math.Log(intensity));
								indices.Add(i);
							}
						}
						double[] ratioSignificanceA = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), 1);
						double[] ratioSignificanceB = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), -1);
						for (int i = 0; i < indices.Count; i++) {
							IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(indices[i]);
							group.AddRatioSignificanceA20(exp, ratioSignificanceA[i]);
							group.AddRatioSignificanceB20(exp, ratioSignificanceB[i]);
						}
					}
				}
				{
					List<double> lratio = new List<double>();
					List<double> lIntensity = new List<double>();
					List<int> indices = new List<int>();
					for (int i = 0; i < combinedData.GetProteinGroupCount(); i++) {
						IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(i);
						double ratio;
						double normalizedRatio;
						int count;
						double dev;
						int nOutliers;
						double pval;
						group.GetRatio21(combinedData, out ratio, out normalizedRatio, out count, out dev, out nOutliers, out pval);
						double intensity = group.GetIntensity0(combinedData) + group.GetIntensity1(combinedData) +
										   group.GetIntensity2(combinedData);
						if (normalizedRatio > 0 && !double.IsInfinity(normalizedRatio) && !double.IsNaN(normalizedRatio) && intensity > 0) {
							lratio.Add(Math.Log(normalizedRatio));
							lIntensity.Add(Math.Log(intensity));
							indices.Add(i);
						}
					}
					double[] ratioSignificanceA = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), 1);
					double[] ratioSignificanceB = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), -1);
					for (int i = 0; i < indices.Count; i++) {
						IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(indices[i]);
						group.RatioSignificanceA21 = ratioSignificanceA[i];
						group.RatioSignificanceB21 = ratioSignificanceB[i];
					}
				}
				if (fileInfo.HasExperiment) {
					string[] exps = fileInfo.GetAllExperimentValues();
					foreach (string exp in exps) {
						List<double> lratio = new List<double>();
						List<double> lIntensity = new List<double>();
						List<int> indices = new List<int>();
						for (int i = 0; i < combinedData.GetProteinGroupCount(); i++) {
							IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(i);
							double ratio;
							double normalizedRatio;
							int count;
							double dev;
							int nOutliers;
							double pval;
							group.GetRatio21(combinedData, out ratio, out normalizedRatio, out count, out dev, out nOutliers,
											 fileInfo.GetRawFileIndicesFromExperiment(exp), out pval);
							double intensity = group.GetIntensity0(combinedData) + group.GetIntensity1(combinedData) +
											   group.GetIntensity2(combinedData);
							if (normalizedRatio > 0 && !double.IsInfinity(normalizedRatio) && !double.IsNaN(normalizedRatio) && intensity > 0) {
								lratio.Add(Math.Log(normalizedRatio));
								lIntensity.Add(Math.Log(intensity));
								indices.Add(i);
							}
						}
						double[] ratioSignificanceA = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), 1);
						double[] ratioSignificanceB = NumUtil.MovingBoxPlot(lratio.ToArray(), lIntensity.ToArray(), -1);
						for (int i = 0; i < indices.Count; i++) {
							IIdentifiedProteinGroup group = combinedData.GetIProteinGroupAt(indices[i]);
							group.AddRatioSignificanceA21(exp, ratioSignificanceA[i]);
							group.AddRatioSignificanceB21(exp, ratioSignificanceB[i]);
						}
					}
				}
			}
		}

		public static void GetRatio10(SilacEvidence[] silacEvidences, bool[] invertSil, IsotopeMsmsEvidence[] isoEvidences,
							  bool[] invertIso, out double ratio,
							  out double normalizedRatio, out int count, out double dev, out int nOutliers,
							  out double pval, int minRatioCount) {
			GetRatio(silacEvidences, invertSil, isoEvidences, invertIso, RatioType.Type10, out ratio, out normalizedRatio,
					 out count, out dev, out nOutliers, out pval, minRatioCount);
		}

		public static void GetRatio20(SilacEvidence[] silacEvidences, bool[] invertSil, IsotopeMsmsEvidence[] isoEvidences,
									  bool[] invertIso, out double ratio,
									  out double normalizedRatio, out int count, out double dev, out int nOutliers,
									  out double pval, int minRatioCount) {
			GetRatio(silacEvidences, invertSil, isoEvidences, invertIso, RatioType.Type20, out ratio, out normalizedRatio,
					 out count, out dev, out nOutliers, out pval, minRatioCount);
		}

		public static void GetRatio21(SilacEvidence[] silacEvidences, bool[] invertSil, IsotopeMsmsEvidence[] isoEvidences,
									  bool[] invertIso, out double ratio,
									  out double normalizedRatio, out int count, out double dev, out int nOutliers,
									  out double pval, int minRatioCount) {
			GetRatio(silacEvidences, invertSil, isoEvidences, invertIso, RatioType.Type21, out ratio, out normalizedRatio,
					 out count, out dev, out nOutliers, out pval, minRatioCount);
		}

		public static void GetRatio(SilacEvidence[] silacEvidence, bool[] invertSil, IsotopeMsmsEvidence[] isoEvidence,
									bool[] invertIso, RatioType type, out double ratio, out double normalizedRatio,
									out int count, out double dev, out int nOutliers, out double pval, int minRatioCount) {
			ratio = 0;
			normalizedRatio = 0;
			count = 0;
			pval = 1;
			List<double> lra = new List<double>();
			List<double> lnra = new List<double>();
			for (int i = 0; i < silacEvidence.Length; i++) {
				SilacEvidence se = silacEvidence[i];
				double r = se.GetRatio(type);
				if (r > 0 && !double.IsInfinity(r) && !double.IsNaN(r)) {
					lra.Add(invertSil[i] ? -Math.Log(r) : Math.Log(r));
					lnra.Add(invertSil[i] ? -Math.Log(se.GetNormalizedRatio(type)) : Math.Log(se.GetNormalizedRatio(type)));
					count++;
				}
			}
			for (int i = 0; i < isoEvidence.Length; i++) {
				IsotopeMsmsEvidence se = isoEvidence[i];
				double r = se.GetRatio(type);
				if (r > 0 && !double.IsInfinity(r) && !double.IsNaN(r)) {
					lra.Add(invertIso[i] ? -Math.Log(r) : Math.Log(r));
					lnra.Add(invertIso[i] ? -Math.Log(se.GetNormalizedRatio(type)) : Math.Log(se.GetNormalizedRatio(type)));
					count++;
				}
			}
			if (lra.Count >= minRatioCount) {
				ratio = ArrayUtil.Median(lra.ToArray());
			}
			if (lra.Count > 1 && lra.Count >= minRatioCount) {
				double dummy;
				new OneSampleTTest().Test(lnra.ToArray(), 0, out dummy, out pval, out dummy, out dummy);
			}
			if (lnra.Count > 0 && lnra.Count >= minRatioCount) {
				normalizedRatio = ArrayUtil.Median((lnra.ToArray()));
			}
			if (count < minRatioCount) {
				ratio = double.NaN;
				normalizedRatio = double.NaN;
				dev = double.NaN;
				nOutliers = 0;
				pval = double.NaN;
				return;
			}
			if (count > 1 && count >= minRatioCount) {
				List<double> lnvals = new List<double>();
				for (int i = 0; i < silacEvidence.Length; i++) {
					SilacEvidence se = silacEvidence[i];
					double nr = se.GetNormalizedRatio(type);
					if (nr > 0 && !double.IsInfinity(nr) && !double.IsNaN(nr)) {
						double lnr = Math.Log(nr);
						lnvals.Add(invertSil[i] ? -lnr : lnr);
					}
				}
				for (int i = 0; i < isoEvidence.Length; i++) {
					IsotopeMsmsEvidence se = isoEvidence[i];
					double nr = se.GetNormalizedRatio(type);
					if (nr > 0 && !double.IsInfinity(nr) && !double.IsNaN(nr)) {
						double lnr = Math.Log(nr);
						lnvals.Add(invertIso[i] ? -lnr : lnr);
					}
				}
				NumUtil.CountOutliers(lnvals.ToArray(), out nOutliers);
				dev = 100 * Math.Sqrt(ArrayUtil.Variance(lnvals.ToArray()));
			} else {
				dev = double.NaN;
				nOutliers = 0;
			}
			ratio = Math.Exp(ratio);
			normalizedRatio = Math.Exp(normalizedRatio);
		}

		public static double GetIntensity0(SilacEvidence[] silacEvidences, IsotopeMsmsEvidence[] isoEvidences) {
			return GetIntensity(silacEvidences, isoEvidences, 0);
		}

		public static double GetIntensity1(SilacEvidence[] silacEvidences, IsotopeMsmsEvidence[] isoEvidences) {
			return GetIntensity(silacEvidences, isoEvidences, 1);
		}

		public static double GetIntensity2(SilacEvidence[] silacEvidences, IsotopeMsmsEvidence[] isoEvidences) {
			return GetIntensity(silacEvidences, isoEvidences, 2);
		}

		public static double GetIntensity(SilacEvidence[] silacEvidences, IsotopeMsmsEvidence[] isoEvidences, int which) {
			double intens = 0;
			foreach (SilacEvidence se in silacEvidences) {
				intens += se.GetIntensity(which);
			}
			foreach (IsotopeMsmsEvidence ie in isoEvidences) {
				if (ie.IsRequantified()) {
					intens += ie.GetIntensity(which);
				} else {
					int silacIndex = ie.GetBestMsmsHit().SilacIndex;
					if (silacIndex == which) {
						intens += ie.Intensity;
					}
				}
			}
			return intens;
		}
	}
}
