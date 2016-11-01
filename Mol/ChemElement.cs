/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using MaxQuant.Num;
using MaxQuant.Util;

namespace MaxQuant.Mol {
	public class ChemElement {
		public static ChemElement[] elements = new ChemElement[]
		                                       	{
		                                       		new ChemElement("H", new double[] {1.0078250321, 2.014101778},
		                                       		                new double[] {99.9885, 0.0115},
		                                       		                1.00794),
		                                       		new ChemElement("Hx", new double[] {1.0078250321, 2.014101778},
		                                       		                new double[] {0.0, 100.0},
		                                       		                2.014101778),
		                                       		new ChemElement("C", new double[] {12.0, 13.0033548378},
		                                       		                new double[] {98.93, 1.07},
		                                       		                12.0107),
		                                       		new ChemElement("Cx", new double[] {12.0, 13.0033548378},
		                                       		                new double[] {0.0, 100.0},
		                                       		                13.0033548378),
		                                       		new ChemElement("N", new double[] {14.0030740052, 15.0001088984},
		                                       		                new double[] {99.632, 0.368},
		                                       		                14.0067),
		                                       		new ChemElement("Nx", new double[] {14.0030740052, 15.0001088984},
		                                       		                new double[] {0.0, 100.0},
		                                       		                15.0001088984),
		                                       		new ChemElement("O", new double[] {15.9949146221, 16.9991315, 17.9991604},
		                                       		                new double[] {99.757, 0.038, 0.205},
		                                       		                15.9994),
		                                       		new ChemElement("Ox", new double[] {15.9949146221, 17.9991604},
		                                       		                new double[] {0.0, 100.0},
		                                       		                17.9991604),
		                                       		new ChemElement("F", new double[] {18.99840320},
		                                       		                new double[] {100.0},
		                                       		                18.9984032),
		                                       		new ChemElement("Na", new double[] {22.98976967},
		                                       		                new double[] {100.0},
		                                       		                22.989770),
		                                       		new ChemElement("P", new double[] {30.97376151},
		                                       		                new double[] {100.0},
		                                       		                30.973761),
		                                       		new ChemElement("S",
		                                       		                new double[]
		                                       		                	{31.97207069, 32.9714585, 33.96786683, 35.96708088},
		                                       		                new double[] {94.93, 0.76, 4.29, 0.02},
		                                       		                32.065),
		                                       		new ChemElement("Sx", new double[] {31.97207069, 33.96786683},
		                                       		                new double[] {0.0, 100.0},
		                                       		                33.96786683),
		                                       		new ChemElement("Cl", new double[] {34.96885271, 36.9659026},
		                                       		                new double[] {75.78, 24.22},
		                                       		                35.453),
		                                       		new ChemElement("K", new double[] {38.9637069, 39.96399867, 40.96182597},
		                                       		                new double[] {93.2581, 0.0117, 6.7302},
		                                       		                39.0983),
		                                       		new ChemElement("Fe",
		                                       		                new double[] {53.9396148, 55.9349421, 56.9353987, 57.9332805},
		                                       		                new double[] {5.845, 91.754, 2.119, 0.282},
		                                       		                55.845),
		                                       		new ChemElement("Cu", new double[] {62.9296011, 64.9277937},
		                                       		                new double[] {69.17, 30.83},
		                                       		                63.546),
		                                       		new ChemElement("Se",
		                                       		                new double[]
		                                       		                	{
		                                       		                		73.9224766, 75.9192141, 76.9199146, 77.9173095, 79.9165218,
		                                       		                		81.9167000
		                                       		                	},
		                                       		                new double[] {0.89, 9.37, 7.63, 23.77, 49.61, 8.73},
		                                       		                78.96),
		                                       		new ChemElement("Br", new double[] {78.9183376, 80.916291},
		                                       		                new double[] {50.69, 49.31},
		                                       		                79.904),
		                                       		new ChemElement("Mo",
		                                       		                new double[]
		                                       		                	{
		                                       		                		91.906810, 93.9050876, 94.9058415, 95.9046789, 96.9060210,
		                                       		                		97.9054078, 99.907477
		                                       		                	},
		                                       		                new double[] {14.84, 9.25, 15.92, 16.68, 9.55, 24.13, 9.63},
		                                       		                95.94),
		                                       		new ChemElement("I", new double[] {126.904468},
		                                       		                new double[] {100.0},
		                                       		                126.90447),
		                                       		new ChemElement("Hg",
		                                       		                new double[]
		                                       		                	{
		                                       		                		195.965815, 197.966752, 198.968262, 199.968309, 200.970285,
		                                       		                		201.970626, 203.973476
		                                       		                	},
		                                       		                new double[] {0.15, 9.97, 16.87, 23.10, 13.18, 29.86, 6.87},
		                                       		                200.59)
		                                       	};

		private readonly double atomicWeight;
		private readonly double[] composition;
		private readonly bool isLabel;
		private readonly double[] masses;
		private readonly double monoisotopicMass;
		private readonly Dictionary<int, double[][]> store = new Dictionary<int, double[][]>();
		private readonly string symbol;

		private ChemElement(string symbol, double[] masses, double[] composition, double atomicWeight) {
			this.symbol = symbol;
			this.masses = masses;
			monoisotopicMass = masses[ArrayUtil.MaxInd(composition)];
			this.composition = composition;
			for (int i = 0; i < this.composition.Length; i++) {
				this.composition[i] *= 0.01;
			}
			this.atomicWeight = atomicWeight;
			isLabel = symbol.EndsWith("x");
		}

		public string Symbol {
			get { return symbol; }
		}

		public bool IsLabel {
			get { return isLabel; }
		}

		public double AtomicWeight {
			get { return atomicWeight; }
		}

		public double MonoIsotopicMass {
			get { return monoisotopicMass; }
		}

		public double[][] GetIsotopeDistribution(int n) {
			if (store.ContainsKey(n)) {
				return store[n];
			}
			double[][] dist = GetIsotopeDistribution(n, masses, composition);
			if (n <= 50) {
				store.Add(n, dist);
			}
			return dist;
		}

		public static double[][] GetIsotopeDistribution(int n, double[] masses, double[] composition) {
			int len = masses.Length;
			int[][] partitions = NumUtil.GetPartitions(n, len);
			double[] ms = new double[partitions.Length];
			double[] weights = new double[partitions.Length];
			for (int i = 0; i < partitions.Length; i++) {
				weights[i] = 1;
				int[] partition = partitions[i];
				for (int j = 0; j < len; j++) {
					ms[i] += partition[j] * masses[j];
					for (int k = 0; k < partition[j]; k++) {
						weights[i] *= composition[j];
					}
				}
				weights[i] *= NumUtil.Multinomial(n, partition);
			}
			int[] o = ArrayUtil.Order(ms);
			ms = ArrayUtil.SubArray(ms, o);
			weights = ArrayUtil.SubArray(weights, o);
			double[][] x = FilterWeights(ms, weights, null, 1e-6);
			ms = x[0];
			weights = x[1];
			x = FilterMasses(ms, weights, null, 0.2);
			ms = x[0];
			weights = x[1];
			return new double[][] {ms, weights};
		}

		public static double[][] FilterMasses(double[] masses, double[] weights, double[] errors, double massPrec) {
			int count = 0;
			double[] newMasses = new double[masses.Length];
			double[] newWeights = new double[masses.Length];
			double[] newErrors = null;
			if (errors != null) {
				newErrors = new double[masses.Length];
			}
			for (int i = 0; i < masses.Length; i++) {
				if (i == masses.Length - 1 ||
				    masses[i + 1] - masses[i] >= massPrec) {
					newMasses[count] = masses[i];
					newWeights[count] = weights[i];
					if (errors != null) {
						newErrors[count] = errors[i];
					}
					count++;
				} else {
					int start = i;
					while (i < masses.Length - 1 &&
					       masses[i + 1] - masses[i] < massPrec) {
						i++;
					}
					for (int j = start; j <= i; j++) {
						newWeights[count] += weights[j];
						newMasses[count] += weights[j] * masses[j];
						if (errors != null) {
							newErrors[count] += weights[j] * errors[j] * errors[j];
						}
					}
					newMasses[count] /= newWeights[count];
					if (errors != null) {
						newErrors[count] /= newWeights[count];
						newErrors[count] = Math.Sqrt(newErrors[count]);
					}
					count++;
				}
			}
			masses = ArrayUtil.SubArray(newMasses, count);
			weights = ArrayUtil.SubArray(newWeights, count);
			if (errors != null) {
				errors = ArrayUtil.SubArray(errors, count);
				return new double[][] {masses, weights, errors};
			}
			return new double[][] {masses, weights};
		}

		public static double[][] FilterWeights(double[] masses, double[] weights, double[] errors, double weightCut) {
			int[] valids = new int[weights.Length];
			int count = 0;
			for (int i = 0; i < valids.Length; i++) {
				if (weights[i] >= weightCut) {
					valids[count++] = i;
				}
			}
			valids = ArrayUtil.SubArray(valids, count);
			if (errors == null) {
				return new double[][]
				       	{
				       		ArrayUtil.SubArray(masses, valids),
				       		ArrayUtil.SubArray(weights, valids)
				       	};
			}
			return new double[][]
			       	{
			       		ArrayUtil.SubArray(masses, valids),
			       		ArrayUtil.SubArray(weights, valids),
			       		ArrayUtil.SubArray(errors, valids)
			       	};
		}

		public int GetMonoIsotopicNucleonCount() {
			return (int) Math.Round(masses[0]);
		}
	}
}