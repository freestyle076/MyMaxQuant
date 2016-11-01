/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.Text;

namespace MaxQuant.Mol {
	public static class MolUtil {
		public static readonly double MassProton = 1.0072764666;

		public static readonly double MassC = CalcMonoMass("C");
		public static readonly double MassH = CalcMonoMass("H");
		public static readonly double MassN = CalcMonoMass("N");
		public static readonly double MassO = CalcMonoMass("O");
		public static readonly double MassP = CalcMonoMass("P");
		public static readonly double MassS = CalcMonoMass("S");
		public static readonly double S34S32Diff = CalcMonoMass("Sx") - CalcMonoMass("S");
		public static readonly double C13C12Diff = CalcMonoMass("Cx") - MassC;
		public static readonly double SulphurShift = 2 * C13C12Diff - S34S32Diff;
		public static readonly double MassNormalCTerminus = MassO + MassH;
		public static readonly double MassNormalNTerminus = MassH;

		public static readonly double WeightC = CalcWeight("C");
		public static readonly double WeightH = CalcWeight("H");
		public static readonly double WeightN = CalcWeight("N");
		public static readonly double WeightO = CalcWeight("O");
		public static readonly double WeightP = CalcWeight("P");
		public static readonly double WeightS = CalcWeight("S");
		public static readonly double WeightNormalCTerminus = WeightO + WeightH;
		public static readonly double WeightNormalNTerminus = WeightH;

		public static readonly double BIonMassOffset = MassProton;
		public static readonly double AIonMassOffset = BIonMassOffset - MassC - MassO;
		public static readonly double YIonMassOffset = MassNormalNTerminus + MassNormalCTerminus + MassProton;
		public static readonly double XIonMassOffset = YIonMassOffset + MassC + MassO;

		public static readonly double BIonWeightOffset = MassProton;
		public static readonly double AIonWeightOffset = BIonWeightOffset - WeightC - WeightO;
		public static readonly double YIonWeightOffset = WeightNormalNTerminus + WeightNormalCTerminus + MassProton;
		public static readonly double XIonWeightOffset = YIonWeightOffset + WeightC + WeightO;

		public static readonly double MassWater = CalcMonoMass("H2O");
		public static readonly double MassAmmonia = CalcMonoMass("NH3");

		private const double AveragineCompositionC = 4.9384;
		private const double AveragineCompositionH = 7.7583;
		private const double AveragineCompositionN = 1.3577;
		private const double AveragineCompositionO = 1.4773;
		private const double AveragineCompositionS = 0.0417;
		private const double binsize = 50;

		private static readonly double AveragineTotal =
			AveragineCompositionC * MassC +
			AveragineCompositionH * MassH +
			AveragineCompositionN * MassN +
			AveragineCompositionO * MassO +
			AveragineCompositionS * MassS;

		private static readonly Dictionary<double, double[][]> isotopePatternListNoSulphur =
			new Dictionary<double, double[][]>();

		private static readonly Dictionary<double, double[][]> isotopePatternListSulphur =
			new Dictionary<double, double[][]>();


		public static readonly double IsotopePatternDiff = GetAverageDifferenceToMonoisotope(1500, 1);

		public static double CalcMonoMass(string formula) {
			return new Molecule(formula).GetMonoIsotopicMass();
		}

		public static double CalcWeight(string formula) {
			return new Molecule(formula).GetMolecularWeight();
		}

		public static double GetAverageDifferenceToMonoisotope(double mass, int peakIndex) {
			if (peakIndex == 0) {
				return 0;
			}
			double[] m = GetAverageIsotopePattern(mass, false)[0];
			if (peakIndex >= m.Length || peakIndex < 0) {
				return peakIndex * (m[1] - m[0]);
			}
			return m[peakIndex] - m[0];
		}

		public static double[][] GetAverageIsotopePattern(double mass, bool includeSulphur) {
			double rounded = Math.Round(mass / binsize) * binsize;
			rounded = Math.Max(rounded, 100);

			if (includeSulphur) {
				lock (isotopePatternListSulphur) {
					if (!isotopePatternListSulphur.ContainsKey(rounded)) {
						isotopePatternListSulphur.Add(rounded,
						                              GetAveraginePeptideMolecule(rounded, includeSulphur).GetIsotopeDistribution(0.2));
					}
					return isotopePatternListSulphur[rounded];
				}
			}
			lock (isotopePatternListNoSulphur) {
				if (!isotopePatternListNoSulphur.ContainsKey(rounded)) {
					isotopePatternListNoSulphur.Add(rounded,
					                                GetAveraginePeptideMolecule(rounded, includeSulphur).GetIsotopeDistribution(0.2));
				}
				return isotopePatternListNoSulphur[rounded];
			}
		}

		public static Molecule GetAveraginePeptideMolecule(double monoMass, bool includeSulphur) {
			return new Molecule(GetAveraginePeptideFormula(monoMass, includeSulphur));
		}

		public static string GetAveraginePeptideFormula(double monoMass, bool includeSulphur) {
			double x;
			if (includeSulphur) {
				x = monoMass / AveragineTotal;
			} else {
				x = monoMass / (AveragineTotal - AveragineCompositionS * MassS);
			}
			int nC = (int) Math.Round(AveragineCompositionC * x);
			int nH = (int) Math.Round(AveragineCompositionH * x);
			int nN = (int) Math.Round(AveragineCompositionN * x);
			int nO = (int) Math.Round(AveragineCompositionO * x);
			int nS = (int) Math.Round(AveragineCompositionS * x);
			StringBuilder result = new StringBuilder();
			if (nC > 0) {
				result.Append("C" + nC);
			}
			if (nH > 0) {
				result.Append("H" + nH);
			}
			if (nN > 0) {
				result.Append("N" + nN);
			}
			if (nO > 0) {
				result.Append("O" + nO);
			}
			if (nS > 0 && includeSulphur) {
				result.Append("S" + nS);
			}
			return result.ToString();
		}
	}
}
