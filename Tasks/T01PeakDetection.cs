/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.IO;
using System.Threading;
using MaxQuant.Num;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Tasks {
	public static class T01PeakDetection {
		public static readonly bool maxIntensity = true;

		public static void Detect(string dir, int missingScans, int pointsForCentroid, CentroidPosition centroidPosition,
		                          bool subtractBackground, int backgroundQuantile, double matchPpm, int minPeaks,
		                          double valleyFactor, bool slicePeaks, double intensityThreshold, IRawFile rawFile,
		                          string basePath, bool maxIntensity, IPeakCollector collector, out double[] centerMassArray,
		                          out float[] centerMassErrorArray, out float[] intensityArray, out float[] minTimeArray,
								  out float[] maxTimeArray, out long[] filePosArray) {
			int npoints = pointsForCentroid;
			FileUtil.Rmdir(dir);
			Directory.CreateDirectory(dir);
			BinaryWriter writer;
			string peaksPath = basePath + ".peaks";
			try {
				writer = FileUtil.GetBinaryWriter(peaksPath);
			} catch (Exception) {
				throw new Exception("Cannot open file " + peaksPath + ". It may be used by another program.");
			}
			List<double> centerMasses = new List<double>();
			List<float> centerMassErrors = new List<float>();
			List<float> intensities = new List<float>();
			List<float> minTimes = new List<float>();
			List<float> maxTimes = new List<float>();
			List<long> filePos = new List<long>();
			int oldMinInd = -1;
			int oldMaxInd = -1;
			List<Ms1CentroidList> cache = new List<Ms1CentroidList>();
			int maxMissingScans = (rawFile.NumberOfMS1MassRanges - 1) * (missingScans + 1) + missingScans;
			double[][] massRanges = new double[rawFile.MS1Count][];
			for (int i = 0; i < rawFile.MS1Count; i++) {
				massRanges[i] = rawFile.GetMS1MassRange(i);
			}
			double[] intensityNorm = CalcIntensityNormalization(rawFile, pointsForCentroid, centroidPosition, subtractBackground,
			                                                    backgroundQuantile, maxIntensity);
			for (int i = 0; i < rawFile.MS1Count; i++) {
				int minInd = Math.Max(0, i - maxMissingScans - 1);
				int maxInd = Math.Min(rawFile.MS1Count - 1, i + maxMissingScans + 1);
				if (i == 0) {
					for (int j = 0; j <= maxInd; j++) {
						Spectrum s = rawFile.GetMS1Spectrum(j, subtractBackground, backgroundQuantile);
						cache.Add(DetectPeaks(s, intensityThreshold, maxIntensity, npoints, centroidPosition));
						s.Dispose();
					}
				} else {
					for (int j = oldMinInd; j < minInd; j++) {
						cache[0].Dispose();
						cache.RemoveAt(0);
					}
					for (int j = oldMaxInd + 1; j <= maxInd; j++) {
						Spectrum s = rawFile.GetMS1Spectrum(j, subtractBackground, backgroundQuantile);
						cache.Add(DetectPeaks(s, intensityThreshold, maxIntensity, npoints, centroidPosition));
						s.Dispose();
					}
				}
				Ms1CentroidList p = cache[i - minInd];
				int[] valids = new int[p.Count];
				int count = 0;
				for (int j = 0; j < p.Count; j++) {
					double cm = p.CenterMass(j);
					bool match = false;
					int ntries = 0;
					for (int k = i - 1; k >= minInd; k--) {
						double[] massRange = massRanges[k];
						if (cm >= massRange[0] && cm <= massRange[1]) {
							Ms1CentroidList q = cache[k - minInd];
							int ind = q.GetClosestIndex(cm);
							double m = q.CenterMass(ind);
							if (ind != -1 && MassMatch(cm, m, matchPpm)) {
								match = true;
								break;
							}
							ntries++;
							if (ntries > missingScans) {
								break;
							}
						}
					}
					if (!match) {
						ntries = 0;
						for (int k = i + 1; k <= maxInd; k++) {
							double[] massRange = massRanges[k];
							if (cm >= massRange[0] && cm <= massRange[1]) {
								Ms1CentroidList q = cache[k - minInd];
								int ind = q.GetClosestIndex(cm);
								double m = q.CenterMass(ind);
								if (ind != -1 && MassMatch(cm, m, matchPpm)) {
									match = true;
									break;
								}
								ntries++;
								if (ntries > missingScans) {
									break;
								}
							}
						}
					}
					if (match) {
						valids[count++] = j;
					}
				}
				valids = ArrayUtil.SubArray(valids, count);
				Ms1CentroidList reduced = p.Extract(valids);
				cache[i - minInd] = reduced;
				byte range = rawFile.GetMS1MassRangeIndex(i);
				for (int j = 0; j < reduced.Count; j++) {
					double cm = reduced.CenterMass(j);
					bool match = false;
					int ntries = 0;
					for (int k = i - 1; k >= minInd; k--) {
						double[] massRange = massRanges[k];
						if (cm >= massRange[0] && cm <= massRange[1]) {
							Ms1CentroidList q = cache[k - minInd];
							int ind = q.GetClosestIndex(cm);
							double m = q.CenterMass(ind);
							if (ind != -1 && MassMatch(cm, m, matchPpm)) {
								GrowablePeak peak = q.GetPeak(ind);
								peak.Add(i, j, reduced, range, intensityNorm[range]);
								match = true;
								break;
							}
							ntries++;
							if (ntries > missingScans) {
								break;
							}
						}
					}
					if (!match) {
						GrowablePeak peak = new GrowablePeak();
						peak.Add(i, j, reduced, range, intensityNorm[range]);
					}
				}
				Ms1CentroidList last = cache[0];
				int nextMinInd = Math.Max(0, i - maxMissingScans);
				for (int j = 0; j < last.Count; j++) {
					GrowablePeak peak = last.GetPeak(j);
					if (peak.IsDisposed()) {
						continue;
					}
					if (peak.LastScanIndex < nextMinInd) {
						Process(peak, writer, centerMasses, centerMassErrors, intensities, minTimes, maxTimes, filePos, minPeaks,
						        rawFile, valleyFactor, slicePeaks, collector, maxIntensity);
						peak.Dispose();
					}
				}
				oldMinInd = minInd;
				oldMaxInd = maxInd;
				if (i % 100 == 0) {
					Thread.Sleep(100);
				}
			}
			for (int i = rawFile.MS1Count - maxMissingScans - 1; i < rawFile.MS1Count; i++) {
				Ms1CentroidList last = cache[i - oldMinInd];
				for (int j = 0; j < last.Count; j++) {
					GrowablePeak peak = last.GetPeak(j);
					if (peak.IsDisposed()) {
						continue;
					}
					if (peak.LastScanIndex == i) {
						Process(peak, writer, centerMasses, centerMassErrors, intensities, minTimes,
						        maxTimes, filePos, minPeaks, rawFile, valleyFactor, slicePeaks, collector, maxIntensity);
						peak.Dispose();
					}
				}
			}
			writer.Close();
			centerMassArray = centerMasses.ToArray();
			int[] o = ArrayUtil.Order(centerMassArray);
			centerMassArray = ArrayUtil.SubArray(centerMassArray, o);
			centerMassErrorArray = ArrayUtil.SubArray(centerMassErrors, o);
			intensityArray = ArrayUtil.SubArray(intensities, o);
			minTimeArray = ArrayUtil.SubArray(minTimes, o);
			maxTimeArray = ArrayUtil.SubArray(maxTimes, o);
			filePosArray = ArrayUtil.SubArray(filePos, o);
		}

		private static bool MassMatch(double m1, double m2, double ppm) {
			double match = 2 * Math.Abs(m1 - m2) / (m1 + m2) * 1e+6;
			return match <= ppm;
		}

		private static Ms1CentroidList DetectPeaks(Spectrum spectrum, double intensityThreshold, bool maxIntensity, int npeaks,
		                                           CentroidPosition centroidPosition) {
			double[] peakCenterMass;
			float[] peakMinMass;
			float[] peakMaxMass;
			float[] peakIntensity;
			DetectPeaks(spectrum, intensityThreshold, out peakCenterMass, out peakMinMass, out peakMaxMass,
			            out peakIntensity, maxIntensity, npeaks, centroidPosition);
			return new Ms1CentroidList(peakCenterMass, peakMinMass, peakMaxMass, peakIntensity);
		}

		public static void DetectPeaks(Spectrum s, double intensityThreshold, out double[] peakCenterMass,
		                               out float[] peakMinMass, out float[] peakMaxMass, out float[] peakIntensity,
		                               bool maxIntensity, int npoints, CentroidPosition centroidPosition) {
			peakCenterMass = new double[s.Count];
			peakMinMass = new float[s.Count];
			peakMaxMass = new float[s.Count];
			peakIntensity = new float[s.Count];
			int peakCount = 0;
			for (int i = 2; i < s.Count - 2; i++) {
				float m2 = s.GetIntensity(i - 2);
				float m1 = s.GetIntensity(i - 1);
				float x = s.GetIntensity(i);
				float p1 = s.GetIntensity(i + 1);
				float p2 = s.GetIntensity(i + 2);
				if (x >= intensityThreshold) {
					if (s.IsMax(x, m1, p1, m2, p2)) {
						int minInd = s.CalcMinPeakIndex(i);
						int maxInd = s.CalcMaxPeakIndex(i);
						if (maxInd - minInd > 2) {
							if (maxInd > i && minInd < i) {
								maxInd--;
								minInd++;
							} else if (maxInd > i) {
								maxInd = i + 1;
							} else if (minInd < i) {
								minInd = i - 1;
							}
							CalcCenterMass(minInd, i, maxInd, s, out peakIntensity[peakCount],
							               out peakCenterMass[peakCount], maxIntensity, npoints, centroidPosition);
							if ((!double.IsNaN(peakCenterMass[peakCount])) && (!float.IsNaN(peakIntensity[peakCount]))) {
								if (minInd > 0) {
									peakMinMass[peakCount] = (float) (0.5 * (s.GetMass(minInd) + s.GetMass(minInd - 1)));
								} else {
									peakMinMass[peakCount] = (float) (1.5 * s.GetMass(0) - 0.5 * s.GetMass(1));
								}
								if (maxInd < s.Count - 1) {
									peakMaxMass[peakCount] = (float) (0.5 * (s.GetMass(maxInd) + s.GetMass(maxInd + 1)));
								} else {
									peakMaxMass[peakCount] = (float) (1.5 * s.GetMass(maxInd) - 0.5 * s.GetMass(maxInd - 1));
								}
								peakCount++;
							}
						}
					}
				}
			}
			peakCenterMass = ArrayUtil.SubArray(peakCenterMass, peakCount);
			peakMinMass = ArrayUtil.SubArray(peakMinMass, peakCount);
			peakMaxMass = ArrayUtil.SubArray(peakMaxMass, peakCount);
			peakIntensity = ArrayUtil.SubArray(peakIntensity, peakCount);
		}

		private static void CalcCenterMass(int minInd, int centerInd, int maxInd, Spectrum s, out float peakIntensity,
		                                   out double peakCenterMass, bool maxIntensity, int npoints,
		                                   CentroidPosition centroidPosition) {
			peakIntensity = 0;
			for (int j = minInd; j <= maxInd; j++) {
				float intensity = s.GetIntensity(j);
				if (maxIntensity) {
					if (intensity > peakIntensity) {
						peakIntensity = intensity;
					}
				} else {
					peakIntensity += intensity;
				}
			}
			if (minInd == maxInd) {
				peakCenterMass = s.GetMass(maxInd);
				return;
			}
			if (minInd == centerInd) {
				peakCenterMass = Estimate2(s.GetMass(centerInd), s.GetMass(centerInd + 1), s.GetIntensity(centerInd),
				                           s.GetIntensity(centerInd + 1));
				return;
			}
			if (maxInd == centerInd) {
				peakCenterMass = Estimate2(s.GetMass(centerInd - 1), s.GetMass(centerInd), s.GetIntensity(centerInd - 1),
				                           s.GetIntensity(centerInd));
				return;
			}
			if (npoints <= 3) {
				switch (centroidPosition) {
					case CentroidPosition.gaussian:
						peakCenterMass = Estimate3(s.GetMass(centerInd - 1), s.GetMass(centerInd), s.GetMass(centerInd + 1),
						                           s.GetIntensity(centerInd - 1), s.GetIntensity(centerInd), s.GetIntensity(centerInd + 1));
						break;
					case CentroidPosition.weightedMean:
						peakCenterMass =
							EstimateWeightedMean(new double[] {s.GetMass(centerInd - 1), s.GetMass(centerInd), s.GetMass(centerInd + 1)},
							                     new double[]
							                     	{s.GetIntensity(centerInd - 1), s.GetIntensity(centerInd), s.GetIntensity(centerInd + 1)});
						break;
					default:
						throw new Exception("Never get here.");
				}
				return;
			}
			int nleft;
			int nright;
			if (npoints % 2 == 1) {
				int d = npoints / 2;
				nleft = Math.Max(centerInd - d, minInd);
				nright = Math.Min(centerInd + d, maxInd);
			} else {
				int d = npoints / 2 - 1;
				nleft = Math.Max(centerInd - d, minInd);
				nright = Math.Min(centerInd + d, maxInd);
				if (nleft != minInd && nright != maxInd) {
					if (s.GetIntensity(nleft - 1) > s.GetIntensity(nright + 1)) {
						nleft--;
					} else {
						nright++;
					}
				} else if (nleft != minInd) {
					nleft--;
				} else if (nright != maxInd) {
					nright++;
				}
			}
			double[] x = new double[nright - nleft + 1];
			double[] y = new double[nright - nleft + 1];
			for (int i = 0; i < x.Length; i++) {
				x[i] = s.GetMass(nleft + i);
				y[i] = s.GetIntensity(nleft + i);
			}
			peakCenterMass = EstimateN(x, y, centerInd - nleft);
			return;
		}

		private static double EstimateWeightedMean(double[] x, double[] y) {
			double m = 0;
			double w = 0;
			for (int i = 0; i < x.Length; i++) {
				m += x[i] * y[i];
				w += y[i];
			}
			return m / w;
		}

		private static double EstimateN(double[] x, double[] y, int ic) {
			double dm = x[ic];
			for (int i = 0; i < y.Length; i++) {
				x[i] -= dm;
				y[i] = Math.Log(y[i]);
			}
			//int im = ic - 1;
			//int ip = ic + 1;
			//double a2 = (y[im] * (x[ip] - x[ic]) + y[ic] * (x[im] - x[ip]) + y[ip] * (x[ic] - x[im])) /
			//            (x[im] * x[im] * (x[ip] - x[ic]) + x[ic] * x[ic] * (x[im] - x[ip]) + x[ip] * x[ip] * (x[ic] - x[im]));
			//double a1 = (y[ic] - y[im] - a2 * (x[ic] * x[ic] - x[im] * x[im])) / (x[ic] - x[im]);
			//double a0 = y[im] - a1 * x[im] - a2 * x[im] * x[im];
			//double[] a = new double[] {a0, a1, a2};
			double[] a = new double[3];
			NumUtil.LinFit2(x, y, a, Qwert);
			return dm - a[1] / a[2] * 0.5;
		}

		public static void Qwert(double x, double[] a) {
			a[0] = 1;
			a[1] = x;
			a[2] = x * x;
		}

		private static double Estimate3(double m1, double m2, double m3, double i1, double i2, double i3) {
			double l1 = Math.Log(i1);
			double l2 = Math.Log(i2);
			double l3 = Math.Log(i3);
			return 0.5 * ((l1 - l2) * (m3 * m3 - m1 * m1) - (l1 - l3) * (m2 * m2 - m1 * m1)) /
			       ((l1 - l2) * (m3 - m1) - (l1 - l3) * (m2 - m1));
		}

		private static double Estimate2(double m1, double m2, double i1, double i2) {
			return (m1 * i1 + m2 * i2) / (i1 + i2);
		}

		private static void Process(GrowablePeak gpeak, BinaryWriter writer, ICollection<double> centerMasses,
		                            ICollection<float> centerMassErrs, ICollection<float> intensities,
		                            ICollection<float> minTimes, ICollection<float> maxTimes,
		                            ICollection<long> filePos, int minPeaks, IRawFile rawFile, double valleyFactor,
									bool split, IPeakCollector collector, bool maxIntensity) {
			Peak peak = gpeak.ToPeak();
			if (peak.Count >= minPeaks) {
				peak.RemoveDoublePeaks(rawFile, maxIntensity);
				if (peak.Count >= minPeaks) {
					peak.Smooth(maxIntensity);
					if (collector != null) {
						collector.AddPeak(peak);
					}
					Peak[] peaks = peak.Decompose(valleyFactor, split, maxIntensity);
					foreach (Peak p in peaks) {
						if (p.Count >= minPeaks) {
							centerMasses.Add(p.CenterMass);
							centerMassErrs.Add(p.CenterMassError);
							intensities.Add(p.Intensity);
							minTimes.Add(p.GetMinTime(rawFile));
							maxTimes.Add(p.GetMaxTime(rawFile));
							filePos.Add(writer.BaseStream.Position);
							p.Write(writer);
						}
						p.Dispose();
					}
				}
			}
		}

		public static void DetectPeaks(Spectrum s, bool maxIntensity, int npeaks, CentroidPosition centroidPosition,
		                               out double[] peakCenterMass, out float[] peakSummedIntensity) {
			float[] peakMinMass;
			float[] peakMaxMass;
			DetectPeaks(s, 0, out peakCenterMass, out peakMinMass, out peakMaxMass, out peakSummedIntensity, maxIntensity,
			            npeaks, centroidPosition);
		}

		private static double[] CalcIntensityNormalization(IRawFile rawFile, int pointsForCentroid,
		                                                   CentroidPosition centroidPosition,
		                                                   bool subtractBackground, int backgroundQuantile, bool maxIntensity) {
			int n = rawFile.NumberOfMS1MassRanges;
			if (n == 0) {
				return new double[] {};
			}
			if (n == 1) {
				return new double[] {1};
			}
			int npoints = pointsForCentroid;
			double[] logAvs = new double[n];
			long[] counts = new long[n];
			for (int i = 0; i < rawFile.MS1Count; i++) {
				byte range = rawFile.GetMS1MassRangeIndex(i);
				Spectrum s = rawFile.GetMS1Spectrum(i, subtractBackground, backgroundQuantile);
				double[] specMasses;
				float[] specIntensities;
				DetectPeaks(s, maxIntensity, npoints, centroidPosition, out specMasses, out specIntensities);
				for (int j = 0; j < specMasses.Length; j++) {
					logAvs[range] += Math.Log(specIntensities[j]);
					counts[range]++;
				}
			}
			for (int i = 0; i < n; i++) {
				logAvs[i] /= counts[i];
			}
			double[] norm = new double[n];
			norm[0] = 1;
			for (int i = 1; i < n; i++) {
				norm[i] = Math.Exp(logAvs[0] - logAvs[i]);
			}
			return norm;
		}
	}
}