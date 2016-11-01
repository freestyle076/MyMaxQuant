/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using MaxQuant.Util;

namespace MaxQuant.Mol {
	public class Molecule {
		public static Molecule Water = new Molecule("H2O");

		private readonly int[] atomCount;
		private readonly int[] atomType;
		private readonly string name;
		private double molecularWeight;

		private double monoIsotopicMass;

		private double mostLikelyMass = Double.NaN;

		public Molecule(string empiricalFormula, string name)
			: this(ProcessEmpiricalFormulaComplex(empiricalFormula), name) {
		}

		public Molecule(int[][] atomTypeAndCount, string name) {
			this.name = name;
			atomType = atomTypeAndCount[0];
			atomCount = atomTypeAndCount[1];
			CalcMasses();
		}

		public Molecule(string empiricalFormula)
			: this(empiricalFormula, empiricalFormula) {
		}

		public int[] AtomType {
			get { return atomType; }
		}

		public int[] AtomCount {
			get { return atomCount; }
		}

		public double GetMonoIsotopicMass() {
			return monoIsotopicMass;
		}

		public double GetMolecularWeight() {
			return molecularWeight;
		}

		public double GetMostLikelyMass(double massPrecision) {
			if (Double.IsNaN(mostLikelyMass)) {
				double[][] distrib = GetIsotopeDistribution(massPrecision);
				double[] masses = distrib[0];
				double[] weights = distrib[1];
				double max = 0;
				int maxind = -1;
				for (int i = 0; i < weights.Length; i++) {
					if (weights[i] > max) {
						max = weights[i];
						maxind = i;
					}
				}
				mostLikelyMass = masses[maxind];
			}
			return mostLikelyMass;
		}

		public string GetName() {
			return name;
		}

		private void CalcMasses() {
			monoIsotopicMass = 0;
			molecularWeight = 0;
			for (int i = 0; i < atomType.Length; i++) {
				int count = atomCount[i];
				int type = atomType[i];
				double mim = ChemElement.elements[type].MonoIsotopicMass;
				double aw = ChemElement.elements[type].AtomicWeight;
				monoIsotopicMass += count * mim;
				molecularWeight += count * aw;
			}
		}

		private static int[][] ProcessEmpiricalFormulaComplex(string formula) {
			if (formula.IndexOf('.') == -1) {
				return ProcessEmpiricalFormula(formula);
			}
			string[] w = formula.Split('.');
			int[] counts = new int[ChemElement.elements.Length];
			for (int i = 0; i < w.Length; i++) {
				int[][] a = ProcessEmpiricalFormula(w[i]);
				int[] type = a[0];
				int[] count = a[1];
				for (int j = 0; j < type.Length; j++) {
					counts[type[j]] += count[j];
				}
			}
			List<int> atomTypes = new List<int>();
			List<int> atomCounts = new List<int>();
			for (int i = 0; i < counts.Length; i++) {
				if (counts[i] > 0) {
					atomTypes.Add(i);
					atomCounts.Add(counts[i]);
				}
			}
			return
				new int[][] {atomTypes.ToArray(), atomCounts.ToArray()};
		}

		private static int[][] ProcessEmpiricalFormula(string formula) {
			formula = StringUtil.RemoveWhitespace(formula);
			if (formula.Length == 0) {
				return new int[][] {new int[0], new int[0]};
			}
			int factor = 1;
			if (formula[0] >= '0' && formula[0] <= '9') {
				int index = 0;
				while (formula[index] >= '0' && formula[index] <= '9') {
					index++;
				}
				factor = Int32.Parse(formula.Substring(0, index));
				formula = formula.Substring(index);
			}
			int[] counts = new int[ChemElement.elements.Length];
			while (formula.Length > 0) {
				formula = ProcessFirst(formula, counts);
			}
			List<int> atomTypes = new List<int>();
			List<int> atomCounts = new List<int>();
			for (int i = 0; i < counts.Length; i++) {
				if (counts[i] > 0) {
					atomTypes.Add(i);
					atomCounts.Add(factor * counts[i]);
				}
			}
			return
				new int[][] {atomTypes.ToArray(), atomCounts.ToArray()};
		}

		private static string ProcessFirst(string formula, int[] counts) {
			for (int i = ChemElement.elements.Length - 1; i >= 0; i--) {
				string symbol = ChemElement.elements[i].Symbol;
				if (formula.ToLower().StartsWith(symbol.ToLower())) {
					formula = formula.Substring(symbol.Length);
					int index = 0;
					while (index < formula.Length && formula[index] >= '0' && formula[index] <= '9') {
						index++;
					}
					int amount;
					if (index == 0) {
						amount = 1;
					} else {
						amount = Int32.Parse(formula.Substring(0, index));
						formula = formula.Substring(index);
					}
					counts[i] += amount;
					return formula;
				}
			}
			throw new Exception("Cannot process " + formula);
		}

		public double[][] GetIsotopeDistribution(double massPrecision) {
			ChemElement element = ChemElement.elements[atomType[0]];
			double[][] distrib = element.GetIsotopeDistribution(atomCount[0]);
			for (int i = 1; i < atomType.Length; i++) {
				element = ChemElement.elements[atomType[i]];
				distrib = Convolute(distrib, element.GetIsotopeDistribution(atomCount[i]), massPrecision, 1e-6);
			}
			return distrib;
		}

		public double[][] GetIsotopeSpectrum(double sigma, int pointsPerSigma, double massPrecision) {
			double spacing = sigma / pointsPerSigma;
			double[][] distrib = GetIsotopeDistribution(massPrecision);
			double[] masses = distrib[0];
			double[] weights = distrib[1];
			double start = masses[0] - 5 * sigma;
			double end = masses[masses.Length - 1] + 5 * sigma;
			double len = end - start;
			int n = (int) Math.Round(len / spacing);
			double[] newMasses = new double[n];
			double[] newWeights = new double[n];
			for (int i = 0; i < n; i++) {
				double mass = start + i * spacing;
				newMasses[i] = mass;
				for (int j = 0; j < masses.Length; j++) {
					newWeights[i] += weights[j] * Math.Exp(-(mass - masses[j]) * (mass - masses[j]) / 2 / sigma / sigma);
				}
			}
			return new double[][] {newMasses, newWeights};
		}

		public static double[][] Convolute(double[][] distrib1, double[][] distrib2,
		                                   double massPrecision, double weightCutoff) {
			double[] masses1 = distrib1[0];
			double[] masses2 = distrib2[0];
			double[] weights1 = distrib1[1];
			double[] weights2 = distrib2[1];
			double[] masses = new double[masses1.Length * masses2.Length];
			double[] weights = new double[masses1.Length * masses2.Length];
			int count = 0;
			for (int i = 0; i < masses1.Length; i++) {
				for (int j = 0; j < masses2.Length; j++) {
					masses[count] = masses1[i] + masses2[j];
					weights[count] = weights1[i] * weights2[j];
					count++;
				}
			}
			int[] o = ArrayUtil.Order(masses);
			masses = ArrayUtil.SubArray(masses, o);
			weights = ArrayUtil.SubArray(weights, o);
			double[][] x = ChemElement.FilterMasses(masses, weights, null, massPrecision);
			masses = x[0];
			weights = x[1];
			if (!Double.IsNaN(weightCutoff)) {
				x = ChemElement.FilterWeights(masses, weights, null, weightCutoff);
				masses = x[0];
				weights = x[1];
			}
			return new double[][] {masses, weights};
		}

		public static void ConvoluteWithErrors(double[] masses1, double[] weights1, double[] errors1, double[] masses2,
		                                       double[] weights2, double massCutoff, double weightCutoff, out double[] masses,
		                                       out double[] weights, out double[] errors) {
			double minMass = masses1[0] + masses2[0];
			double maxMass = masses1[masses1.Length - 1] + masses2[masses2.Length - 1];
			int nbins = (int) Math.Ceiling((maxMass - minMass) / massCutoff) + 1;
			double[] m = new double[nbins];
			double[] w = new double[nbins];
			double[] e = new double[nbins];
			for (int i = 0; i < masses1.Length; i++) {
				double error = errors1[i];
				double e2 = error * error;
				for (int j = 0; j < masses2.Length; j++) {
					double mass = masses1[i] + masses2[j];
					double weight = weights1[i] * weights2[j];
					int bin = (int) Math.Round((mass - minMass) / massCutoff);
					m[bin] += mass * weight;
					w[bin] += weight;
					e[bin] += e2 * weight;
				}
			}
			List<double> mm = new List<double>();
			List<double> ww = new List<double>();
			List<double> ee = new List<double>();
			for (int i = 0; i < nbins; i++) {
				if (m[i] != 0 && w[i] >= weightCutoff) {
					mm.Add(m[i] / w[i]);
					ww.Add(w[i]);
					ee.Add(Math.Sqrt(e[i] / w[i]));
				}
			}
			masses = mm.ToArray();
			weights = ww.ToArray();
			errors = ee.ToArray();
		}

		public static void ConvoluteWithErrorsRecent(double[] masses1, double[] weights1, double[] errors1, double[] masses2,
		                                             double[] weights2,
		                                             double massCutoff, double weightCutoff, out double[] masses,
		                                             out double[] weights, out double[] errors) {
			double minMass = masses1[0] + masses2[0];
			double maxMass = masses1[masses1.Length - 1] + masses2[masses2.Length - 1];
			int nbins = (int) Math.Ceiling((maxMass - minMass) / massCutoff) + 1;
			List<double>[] m = new List<double>[nbins];
			List<double>[] w = new List<double>[nbins];
			List<double>[] e = new List<double>[nbins];
			for (int i = 0; i < masses1.Length; i++) {
				for (int j = 0; j < masses2.Length; j++) {
					double mass = masses1[i] + masses2[j];
					double weight = weights1[i] * weights2[j];
					double error = errors1[i];
					int bin = (int) Math.Round((mass - minMass) / massCutoff);
					if (m[bin] == null) {
						m[bin] = new List<double>();
						w[bin] = new List<double>();
						e[bin] = new List<double>();
					}
					m[bin].Add(mass);
					w[bin].Add(weight);
					e[bin].Add(error);
				}
			}
			List<double> mm = new List<double>();
			List<double> ww = new List<double>();
			List<double> ee = new List<double>();
			for (int i = 0; i < nbins; i++) {
				if (m[i] != null) {
					double mmm = 0;
					double www = 0;
					double eee = 0;
					for (int j = 0; j < m[i].Count; j++) {
						www += w[i][j];
						mmm += w[i][j] * m[i][j];
						eee += w[i][j] * e[i][j] * e[i][j];
					}
					if (www >= weightCutoff) {
						mmm /= www;
						eee /= www;
						eee = Math.Sqrt(eee);
						mm.Add(mmm);
						ww.Add(www);
						ee.Add(eee);
					}
				}
			}
			masses = mm.ToArray();
			weights = ww.ToArray();
			errors = ee.ToArray();
		}

		public static Molecule[] GetDifferences(Molecule molecule1, Molecule molecule2) {
			int[] counts1 = new int[ChemElement.elements.Length];
			for (int j = 0; j < molecule1.atomType.Length; j++) {
				counts1[molecule1.atomType[j]] = molecule1.atomCount[j];
			}
			int[] counts2 = new int[ChemElement.elements.Length];
			for (int j = 0; j < molecule2.atomType.Length; j++) {
				counts2[molecule2.atomType[j]] = molecule2.atomCount[j];
			}
			for (int i = 0; i < counts1.Length; i++) {
				if (counts1[i] > counts2[i]) {
					counts1[i] -= counts2[i];
					counts2[i] = 0;
				} else {
					counts2[i] -= counts1[i];
					counts1[i] = 0;
				}
			}
			int[] types1 = new int[counts1.Length];
			int count1 = 0;
			for (int i = 0; i < counts1.Length; i++) {
				if (counts1[i] > 0) {
					types1[count1++] = i;
				}
			}
			types1 = ArrayUtil.SubArray(types1, count1);
			Molecule diff1 = new Molecule(new int[][] {types1, ArrayUtil.SubArray(counts1, types1)}, "");
			int[] types2 = new int[counts2.Length];
			int count2 = 0;
			for (int i = 0; i < counts2.Length; i++) {
				if (counts2[i] > 0) {
					types2[count2++] = i;
				}
			}
			types2 = ArrayUtil.SubArray(types2, count2);
			Molecule diff2 = new Molecule(new int[][] {types2, ArrayUtil.SubArray(counts2, types2)}, "");
			return new Molecule[] {diff1, diff2};
		}

		public static Molecule Sum(Molecule molecule1, Molecule molecule2) {
			return Sum(new Molecule[] {molecule1, molecule2}, new int[] {1, 1});
		}

		public static Molecule Sum(Molecule[] molecules, int[] n) {
			int[] counts = new int[ChemElement.elements.Length];
			for (int i = 0; i < molecules.Length; i++) {
				Molecule mol = molecules[i];
				for (int j = 0; j < mol.atomType.Length; j++) {
					counts[mol.atomType[j]] += n[i] * mol.atomCount[j];
				}
			}
			int[] types = new int[counts.Length];
			int count = 0;
			for (int i = 0; i < counts.Length; i++) {
				if (counts[i] > 0) {
					types[count++] = i;
				}
			}
			types = ArrayUtil.SubArray(types, count);
			return new Molecule(new int[][] {types, ArrayUtil.SubArray(counts, types)}, "");
		}

		public static Molecule Max(Molecule x, Molecule y) {
			int[] counts1 = new int[ChemElement.elements.Length];
			for (int j = 0; j < x.atomType.Length; j++) {
				counts1[x.atomType[j]] = x.atomCount[j];
			}
			int[] counts2 = new int[ChemElement.elements.Length];
			for (int j = 0; j < y.atomType.Length; j++) {
				counts2[y.atomType[j]] = y.atomCount[j];
			}
			int[] counts = new int[ChemElement.elements.Length];
			for (int i = 0; i < ChemElement.elements.Length; i++) {
				counts[i] = Math.Max(counts1[i], counts2[i]);
			}
			int[] types = new int[counts.Length];
			int count = 0;
			for (int i = 0; i < counts.Length; i++) {
				if (counts[i] > 0) {
					types[count++] = i;
				}
			}
			types = ArrayUtil.SubArray(types, count);
			return new Molecule(new int[][] {types, ArrayUtil.SubArray(counts, types)}, "");
		}
	}
}