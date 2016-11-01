/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using MaxQuant.Mol;
using MaxQuant.Spec;
using MaxQuant.Util;

namespace MaxQuant.Tasks {
	public class T02DeIsotoping {
		private List<byte> resultClusterCharges = new List<byte>();
		private List<int[]> resultClusters = new List<int[]>();
		private float[][] tmpCorr;
		private float[] tmpIntens;
		private float[] tmpMassErrors;
		private double[] tmpMasses;
		private Peak[] tmpPeaks;
		private int[][] tmpx;

		
		public static List<int[]> FindIsotopeClusters(
			byte minCharge,
			byte maxCharge, 
			double correlationThreshold,
			double sigmas, 
			double isotopeValleyFactor, 
			out byte[] clusterCharges,
		    int peakCount, 
		    double[] centerMz, 
		    float[] centerMzErrors,
			float[] minTimes, 
			float[] maxTimes, 
			Peak[] peaks, 
			float[] intensities,
		    IPeakList peakList) 
		{

			//** cluster assignments **
			List<int[]> clusterInd = CalcClusterIndices(minCharge, maxCharge, correlationThreshold, peakCount, centerMz,
			                                            centerMzErrors, minTimes, maxTimes, peaks, peakList);
			
			// ** discard any unclustered peaks ** 
			if (peaks != null) {
				for (int i = 0; i < peaks.Length; i++) {
					if (peaks[i] != null) {
						peaks[i].Dispose();
						peaks[i] = null;
					}
				}
			}


			clusterCharges = new byte[clusterInd.Count];
			

			//** this is some of the worst programming I've ever seen **
			T02DeIsotoping deiso = new T02DeIsotoping();

			//** while(true) is this cool new thing... **
			//** or you know, perhaps while(hasChanged)... **

			//** iteratively improve cluster indices until no more improvements **
			for (;;) 
			{
				bool hasChanged;

				
				clusterInd = deiso.ImproveIsotopeClusterIndices(clusterInd, minCharge, maxCharge, correlationThreshold, sigmas,
				                                                isotopeValleyFactor, ref clusterCharges, out hasChanged, centerMz,
				                                                centerMzErrors, intensities, peakList);
				if (!hasChanged) {
					break;
				}
			}
			deiso.Dispose();
			return clusterInd;
		}

		public List<int[]> ImproveIsotopeClusterIndices(IList<int[]> clusterInd, byte minCharge, byte maxCharge,
		                                                double correlationThreshold, double sigmas, double isotopeValleyFactor,
		                                                ref byte[] clusterCharges, out bool hasChanged, double[] centerMz,
		                                                float[] centerMzErrors, float[] intensities, IPeakList peakList) 
		{
			hasChanged = false;
			resultClusters.Clear();
			resultClusterCharges.Clear();

			//** iterate through the clusters **
			for (int i = 0; i < clusterInd.Count; i++) 
			{
				//** initially, all charges are 0 **
				byte charge = clusterCharges[i];

				//** if the charge has already been assigned **
				if (charge > 0) 
				{
					//** if charge is outside of prescribed range **
					//** how does a cluster get assigned an out-of-bounds charge? **
					if (charge < minCharge || charge > maxCharge) 
					{
						hasChanged = true;
					} 
					else 
					{
						resultClusters.Add(clusterInd[i]);
						resultClusterCharges.Add(charge);
					}
				} 

				//** else the cluster has not been assigned a charge **
				else 
				{
					byte ch;
					//** tmpx is a 2D array of ints **
					//** could be [original cluster] **
					//** or [new cluster],[those removed] **
					tmpx =
						SplitIsotopeCluster(clusterInd[i], minCharge, maxCharge, correlationThreshold, sigmas, isotopeValleyFactor, out ch,
						                    centerMz, centerMzErrors, intensities, peakList);

					//** if return has two arrays [new cluster][those removed] **
					if (tmpx.Length > 1) 
					{
						// it's changed
						hasChanged = true;
					}

					//** iterate through each returned array **
					for (int j = 0; j < tmpx.Length; j++) 
					{
						//** a returned array with 1 peak cannot be processed further **
						if (tmpx[j].Length == 1) 
						{
							//do nothing
						} 

						//** else it has more than one peak **
						else if (tmpx[j].Length > 1) 
						{
							//** collect it and its charge **
							resultClusters.Add(tmpx[j]);
							if (j == 0) 
							{
								resultClusterCharges.Add(ch);
							} 
							else 
							{
								resultClusterCharges.Add(0);
							}
						}
						tmpx[j] = null;
					}
				}
				clusterInd[i] = null;
			}

			//** IGNORE**
			for (int i = 0; i < clusterInd.Count; i++) 
			{
				clusterInd[i] = null;
			}
			clusterInd.Clear();
			clusterCharges = resultClusterCharges.ToArray();
			return new List<int[]>(resultClusters);
		}

		//** splits a cluster into the longest envelope-consistent cluster and the excluded peaks **
		private int[][] SplitIsotopeCluster(int[] m, byte minCharge, byte maxCharge, double threshold, double sigmas,
		                                    double isotopeValleyFactor, out byte resultCharge, double[] centerMz,
		                                    float[] centerMzErrors, float[] intensities, IPeakList peakList) {
			//** !!!!!!!!!! m is cluster !!!!!!!!!  **
			//** wonderful naming conventions **
			int n = m.Length;
			if (n < 2) {
				throw new Exception("Cannot split cluster of length " + n + ".");
			}


			int[] bestMatches = new int[0];

			//** if the cluster length is over 100 (why 100?)**
			if (n > 100) 
			{
				//TODO: find more elegant way of splitting
				resultCharge = 0;
				bestMatches = new int[n / 2];
				for (int i = 0; i < bestMatches.Length; i++) 
				{
					bestMatches[i] = i;
				}
			} 


			else 
			{
				//** cosine correlations between each peak and all left-hand peaks **
				tmpCorr = new float[n][];
				
				//** stores all peaks in cluster locally **
				tmpPeaks = new Peak[n];


				//** get all peaks in cluster **
				for (int i = 0; i < n; i++) 
				{
					tmpPeaks[i] = peakList.GetPeakDiscard(m[i]);
				}

				//** calculate the correlation between the each peak and all left-hand peaks **
				for (int i = 0; i < n; i++) 
				{
					tmpCorr[i] = new float[i];
					for (int j = 0; j < i; j++) 
					{
						tmpCorr[i][j] = CalcCorrelation(tmpPeaks[i], tmpPeaks[j]);
					}
				}

				//** mass, mass error, intensity of each peak in cluster **
				tmpMasses = ArrayUtil.SubArray(centerMz, m);
				tmpMassErrors = ArrayUtil.SubArray(centerMzErrors, m);
				tmpIntens = ArrayUtil.SubArray(intensities, m);


				byte bestMatchesCharge = 0;
				//** iterate through charge states, descending from maximum **
				for (byte ch = maxCharge; ch >= minCharge; ch--) 
				{
					//** iterate through the cluster **
					for (int i = 0; i < n; i++) 
					{
						//** threshold is correlation threshold **
						//** i is the current peak in the cluster **

						//** get set of peaks in cluster consistent with charge state, and envelope shape **   
						int[] matches =
							GetMatches(ch, i, threshold, sigmas, tmpCorr, tmpMasses, tmpMassErrors, tmpIntens,
							           isotopeValleyFactor);
						
						//** keep the longest match and its charge **
						if (matches.Length > bestMatches.Length) 
						{
							bestMatches = matches;
							bestMatchesCharge = ch;
						}
					}
				}

				//** charge of best match **
				resultCharge = bestMatchesCharge;
				
				//** IGNORE **
				for (int i = 0; i < n; i++) 
				{
					tmpCorr[i] = null;
				}
				for (int i = 0; i < tmpPeaks.Length; i++) {
					//tmpPeaks[i].Dispose();
					tmpPeaks[i] = null;
				}
				tmpCorr = null;
				tmpPeaks = null;
				tmpMasses = null;
				tmpMassErrors = null;
				tmpIntens = null;
			}


			bestMatches = ArrayUtil.UniqueValues(bestMatches);

			//** if the best match is shorter than the inputted cluster **
			if (bestMatches.Length < n) 
			{
				//** then return the matched set and those not in the best match **
				return new int[][]
				       	{
				       		ArrayUtil.SubArray(m, bestMatches),
				       		ArrayUtil.SubArray(m, ArrayUtil.Complement(bestMatches, n))
				       	};
			}

			//** else return the original cluster as it is the best match **
			return new int[][] {m};
		}

		//** cluster the peaks **
		public static List<int[]> CalcClusterIndices(int minCharge, int maxCharge, double correlationThreshold, int peakCount,
		                                             double[] centerMz, float[] centerMzErrors, float[] minTimes,
		                                             float[] maxTimes,
		                                             Peak[] peaks, IPeakList peakList) {


			NeighbourList neighbourList = new NeighbourList();

			//** iterate through all peaks **
			for (int j = 0; j < peakCount; j++) 
			{

				//** current peak's mz and RT range **
				double massJ = centerMz[j];
				float massErrorJ = centerMzErrors[j];
				float timeMinJ = minTimes[j];
				float timeMaxJ = maxTimes[j];

				//** get index of nearest peak at or above current mass -1.1 **
				int start = ArrayUtil.CeilIndex(centerMz, massJ - 1.1);
				
				//** get index of nearest peak at or below current mass -1.2 
				int w = ArrayUtil.FloorIndex(centerMz, massJ - 1.2);

				//** remove any peaks outside of massj - 1.2 ** 
				//** so there's a removal to the left of this peak outside 1.2 away... **
				//** what is this all about?? **
				
				for (int i = 0; i < w; i++) {
					if (peaks != null && peaks[i] != null) {
						peaks[i].Dispose();
						peaks[i] = null;
					}
				}

				//** iterate from current peak at mass - 1.1 to current peak **
				//** iterates through left "adjacent" traces, neihbors with current trace (j) if valid **
				for (int i = start; i < j; i++) 
				{

					//** comparing peak mz and RT range **
					double massI = centerMz[i];
					double massErrorI = centerMzErrors[i];
					double timeMinI = minTimes[i];
					double timeMaxI = maxTimes[i];
					
					//** difference in mass and synthesized mass error
					double massDiff = Math.Abs(massI - massJ);
					double massError = 5 * Math.Sqrt(massErrorI * massErrorI + massErrorJ * massErrorJ);
					
					//** invalidating conditions: **
					
					//** 1) mass difference is greater than minimum **
					if (massDiff > MolUtil.C13C12Diff + massError) {
						continue;
					}

					//** 2) no RT overlap
					if (timeMinI >= timeMaxJ) {
						continue;
					}

					//** 2) no RT overlap
					if (timeMinJ >= timeMaxI) {
						continue;
					}

					//** 3) mass difference doesn't match any charge states **
					if (!FitsMassDifference(massDiff, massError, minCharge, maxCharge)) {
						continue;
					}

					//** 4) The intensity profile correlation (cosine similarity) fails the threshold **
					if (CalcCorrelation(peakList.GetPeakKeep(i), peakList.GetPeakKeep(j)) < correlationThreshold) {
						continue;
					}

					//** create an edge between peak I and peak J if valid: **
					//** 1) mass difference exceeds minimum
					//** 2) RT has overlap
					//** 3) mass difference fits a charge state
					//** 4) intensity profiles have strong correlation
					neighbourList.Add(i, j);
				}
			}

			//** convert edge list to clusters! **
			List<int[]> clusterList = new List<int[]>();
			//** iterate through all peaks **
			for (int i = 0; i < peakCount; i++) 
			{
				//** if the peak has neighbors... **
				if (!neighbourList.IsEmptyAt(i)) 
				{
					
					HashSet<int> currentCluster = new HashSet<int>();

					AddNeighbors(i, currentCluster, neighbourList);
					int[] c = SortByMass(currentCluster.ToArray(), centerMz);
					clusterList.Add(c);
				}
			}
			return clusterList;
		}

		private static int[] SortByMass(int[] c, double[] centerMz) {
			double[] m = new double[c.Length];
			for (int j = 0; j < c.Length; j++) {
				m[j] = centerMz[c[j]];
			}
			int[] o = ArrayUtil.Order(m);
			return ArrayUtil.SubArray(c, o);
		}


		//** checks if the mass difference is acceptable for any of the given charge states **
		public static bool FitsMassDifference(double massDiff, double massError, int minCharge, int maxCharge) 
		{
			//** iterate through charge states **
			for (int i = minCharge; i <= maxCharge; i++) 
			{

				double error = Math.Sqrt(
					MolUtil.SulphurShift * MolUtil.SulphurShift / i / i
					+ massError * massError);
				
				if (Math.Abs(massDiff - MolUtil.IsotopePatternDiff / i) <= error) {
					return true;
				}
			}
			return false;
		}

		//** aligns the two peak's intensity profiles, returns cosine similarity (A dot B) / (|A| * |B|) **
		public static float CalcCorrelation(Peak psI, Peak psJ) 
		{
			//** scan indices of each peak **
			int[] indI = psI.GetScanIndices();
			int[] indJ = psJ.GetScanIndices();
			
			//** min and max scan indexes over both peaks **
			int minIndex = Math.Min(indI[0], indJ[0]);
			int maxIndex = Math.Max(indI[indI.Length - 1], indJ[indJ.Length - 1]);
			
			//** range of synthesized indices **
			int len = maxIndex - minIndex + 1;
			
			//** intensity profile of peakI **
			float[] profileI = new float[len + 2];
			for (int a = 0; a < indI.Length; a++) {
				profileI[indI[a] - minIndex + 1] = psI.GetSmoothIntensity(a);
			}
			
			//** " peakJ **
			float[] profileJ = new float[len + 2];
			for (int a = 0; a < indJ.Length; a++) {
				profileJ[indJ[a] - minIndex + 1] = psJ.GetSmoothIntensity(a);
			}

			int c = 0;
			int l = profileI.Length;
			int[] valids = new int[l];

			//** iterate through I's intensity profile, collecting index of valid peak intensities **
			//** valid === one peak has an entry, or its the begnning or end **
			for (int x = 0; x < l; x++) {
				
				if (profileI[x] > 0 || profileJ[x] > 0 || x == 0 || x == l - 1) {
					valids[c++] = x;
				}
			}

			//** trim valid array **
			valids = ArrayUtil.SubArray(valids, c);

			//** return cosine similarity between the two aligned arrays **
			return ArrayUtil.Cosine(ArrayUtil.SubArray(profileI, valids), ArrayUtil.SubArray(profileJ, valids));
		}


		//** i: current peak index **
		//** currentCluster: side effect **
		//** neighbourList: edge set **
		//** creates a cluster  **
		public static void AddNeighbors(int i, HashSet<int> currentCluster,
		                                NeighbourList neighbourList) 
		{

			//** nodes to process **
			Stack<int> todo = new Stack<int>();
			todo.Push(i);


			while (todo.Count > 0) 
			{
				int next = todo.Pop();
				
				//** is the node already in the cluster? **
				if (!currentCluster.Contains(next)) 
				{

					//** collect the current node **
					currentCluster.Add(next);

					//** does the current node have neighbors remaining? **
					if (!neighbourList.IsEmptyAt(next)) 
					{
						foreach (int x in neighbourList.GetAt(next)) 
						{
							todo.Push(x);
						}
						neighbourList.DeleteAt(next);
					}
				}
			}
		}

		public static int[] GetMatches(int charge, int startIndex, double threshold,
		                               double sigmas, float[][] corr, double[] masses,
		                               float[] massErrors, float[] intensities,
		                               double isotopeValleyFactor) 
		{
			//** start index is the index of the current peak **
			//** treshold is the correlation threshold **
			//** corr is the correlation matrix **
			//** sigmas? **
			//** 

			//** mass of the current peak **
			double startMass = masses[startIndex];

			//** right-hand peaks that charge state separation **
			List<int> upMatches = new List<int>();
			
			//** search to the right by charge distance separation **
			//** searching for traces in the envelope **
			for (int i = 1;; i++) 
			{
				//** my god naming!!!!! **

				//** theoretical mass of next trace in envelope **
				double m = startMass + i * MolUtil.IsotopePatternDiff / charge;

				//** find indexes of peaks that fit the description of "m" **
				int[] fits = CollectFittingMasses(startIndex, m, masses, massErrors, sigmas, charge);
				
				//** if no fitting peaks then quit **
				if (fits.Length == 0) 
				{
					break;
				}

				//** this block is unnecessary **
				//** if there's only one fitting mass **
				if (fits.Length == 1) 
				{
					//** get correlation between current peak and matched peak **
					int minInd = Math.Min(startIndex, fits[0]);
					int maxInd = Math.Max(startIndex, fits[0]);
					double cor = corr[maxInd][minInd];
					
					//** if surpasses the threshold collect as match **
					if (cor >= threshold) 
					{
						upMatches.Add(fits[0]);
					} 
					//** else quit searching **
					else 
					{
						break;
					}
				}
				//** else more than one peak fit the description of "m" ** 
				else 
				{
					//** correlations of each fitted peak **
					double[] corrs = new double[fits.Length];
					for (int j = 0; j < fits.Length; j++) 
					{
						int minInd = Math.Min(startIndex, fits[j]);
						int maxInd = Math.Max(startIndex, fits[j]);
						corrs[j] = corr[maxInd][minInd];
					}

					//** find the highest correlated peak at this m/z **
					int[] o = ArrayUtil.Order(corrs);
					int maxIndex = o[o.Length - 1];
					if (corrs[maxIndex] >= threshold) 
					{
						upMatches.Add(fits[maxIndex]);
					} 
					//** if highest correlated fails the threshold, quit searching **
					else 
					{
						break;
					}
				}
			}


			//** REPEAT BUT TO THE LEFT **

			//** left-hand peaks that match current peak **
			List<int> downMatches = new List<int>();
			
			//** iterate through isotopic envelope positions **
			for (int i = 1;; i++) 
			{
				double m = startMass - i * MolUtil.IsotopePatternDiff / charge;
				int[] fits = CollectFittingMasses(startIndex, m, masses, massErrors, sigmas, charge);
				if (fits.Length == 0) {
					break;
				}
				if (fits.Length == 1) {
					int minInd = Math.Min(startIndex, fits[0]);
					int maxInd = Math.Max(startIndex, fits[0]);
					double cor = corr[maxInd][minInd];
					if (cor >= threshold) {
						downMatches.Add(fits[0]);
					} else {
						break;
					}
				} else {
					double[] corrs = new double[fits.Length];
					for (int j = 0; j < fits.Length; j++) {
						int minInd = Math.Min(startIndex, fits[j]);
						int maxInd = Math.Max(startIndex, fits[j]);
						corrs[j] = corr[maxInd][minInd];
					}
					int[] o = ArrayUtil.Order(corrs);
					int maxIndex = o[o.Length - 1];
					if (corrs[maxIndex] >= threshold) {
						downMatches.Add(fits[maxIndex]);
					} else {
						break;
					}
				}
			}

			//** synthesized indexes from left and right-hand matches **
			int[] result = new int[upMatches.Count + downMatches.Count + 1];
			
			//** left-hand **
			for (int i = 0; i < downMatches.Count; i++) 
			{
				result[i] = downMatches[downMatches.Count - i - 1];
			}

			//** current **
			result[downMatches.Count] = startIndex;
			
			//** right-hand **
			for (int i = 0; i < upMatches.Count; i++) 
			{
				result[downMatches.Count + 1 + i] = upMatches[i];
			}

			//** aggregate intensity values of matched peaks **
			float[] profile = ArrayUtil.SubArray(intensities, result);
			
			//** if the current peak is the local minimum return nada **
			if (IsLocalMinimum(profile, downMatches.Count)) 
			{
				return new int[0];
			}

			//** indexes of peaks that are less in intensity than adjacent peaks
			int[] minima = GetLocalMinima(profile, isotopeValleyFactor);
			

			//** if found a local minimum **
			if (minima.Length > 0) 
			{

				//** find lower minimum and upper minimum **
				int lower = 0;
				for (int i = 0; i < minima.Length; i++) 
				{
					if (minima[i] < downMatches.Count) {
						lower = minima[i];
					}
				}
				int upper = result.Length - 1;
				for (int i = minima.Length - 1; i >= 0; i--) {
					if (minima[i] > downMatches.Count) {
						upper = minima[i];
					}
				}
				
				//** prune matching peaks to the upper and lower minima **
				//** envelopes cannot have a valley **
				result = ArrayUtil.SubArray(result, lower, upper + 1);
			}


			//** if pruned down to only 2 and the mass of the current peak > 2500 return current peak only **
			if (charge * startMass > 2500 && result.Length == 2) 
			{
				return new int[] {startIndex};
			}

			//** if mass below 1000 (what's special about mass of 1000?) **
			if (charge * startMass < 1000) 
			{
				//** intensity profile of pruned envelope **
				float[] p = ArrayUtil.SubArray(intensities, result);
				double max = -double.MaxValue;
				int maxPos = -1;

				//** find max intensity **
				for (int i = 0; i < p.Length; i++) 
				{
					if (p[i] > max) 
					{
						max = p[i];
						maxPos = i;
					}
				}

				//** trim results such that most intense is on left **
				if (maxPos != 0) 
				{
					return ArrayUtil.SubArray(result, maxPos, result.Length);
				}
			}
			return result;
		}

		private static int[] CollectFittingMasses(int index, double m, double[] masses,
		                                          float[] massErrors, double sigmas, int charge) {
			List<int> fits = new List<int>();
			for (int i = 0; i < masses.Length; i++) {
				if (i == index) {
					continue;
				}
				double massErrorSq = sigmas * sigmas * (massErrors[index] * massErrors[index]
				                                        + massErrors[i] * massErrors[i]);

				double error = Math.Sqrt(
					MolUtil.SulphurShift * MolUtil.SulphurShift / charge / charge
					+ massErrorSq);
				if (Math.Abs(m - masses[i]) <= error) {
					fits.Add(i);
				}
			}
			return fits.ToArray();
		}

		//** finds indexes of peaks that are lesser in intensity than left and right adjacent peaks **
		private static int[] GetLocalMinima(float[] profile, double isotopeValleyFactor) 
		{
			List<int> mins = new List<int>();
			for (int i = 1; i < profile.Length - 1; i++) 
			{
				if (IsLocalMinimum(profile, i)) 
				{
					mins.Add(i);
				}
			}
			List<int> result = new List<int>();
			foreach (int min in mins) {
				double maxL = GetLeftMax(profile, min);
				double maxR = GetRightMax(profile, min);
				double minval = profile[min];
				double smallMax = Math.Min(maxL, maxR);
				if (smallMax / minval >= isotopeValleyFactor) {
					result.Add(min);
				}
			}
			return result.ToArray();
		}

		private static double GetLeftMax(float[] profile, int min) {
			double max = -Double.MaxValue;
			for (int i = 0; i < min; i++) {
				if (profile[i] > max) {
					max = profile[i];
				}
			}
			return max;
		}

		private static double GetRightMax(float[] profile, int min) {
			double max = -Double.MaxValue;
			for (int i = min + 1; i < profile.Length; i++) {
				if (profile[i] > max) {
					max = profile[i];
				}
			}
			return max;
		}

		private static bool IsLocalMinimum(float[] profile, int index) {
			if (index == 0 || index == profile.Length - 1) {
				return false;
			}
			return profile[index - 1] > profile[index] &&
			       profile[index + 1] > profile[index];
		}

		public void Dispose() {
			if (resultClusterCharges != null) {
				resultClusterCharges.Clear();
				resultClusterCharges = null;
			}
			if (resultClusters != null) {
				resultClusters.Clear();
				resultClusters = null;
			}
			if (tmpCorr != null) {
				for (int i = 0; i < tmpCorr.Length; i++) {
					tmpCorr[i] = null;
				}
				tmpCorr = null;
			}
			if (tmpx != null) {
				for (int i = 0; i < tmpx.Length; i++) {
					tmpx[i] = null;
				}
				tmpx = null;
			}
			tmpIntens = null;
			tmpMassErrors = null;
			tmpMasses = null;
			tmpPeaks = null;
		}
	}
}