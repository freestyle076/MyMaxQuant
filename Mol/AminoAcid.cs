/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using MaxQuant.Util;

namespace MaxQuant.Mol {
	public class AminoAcid : Molecule {
		public static readonly AminoAcid Alanine = new AminoAcid("C3H5NO", "Alanine", "Ala", 'A', 7.4,
																 new string[] { "GCT", "GCC", "GCA", "GCG" });

		public static readonly AminoAcid Arginine = new AminoAcid("C6H12N4O", "Arginine", "Arg", 'R', 4.2,
		                                                          new string[] {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"});

		public static readonly AminoAcid Asparagine = new AminoAcid("C4H6N2O2", "Asparagine", "Asn", 'N', 4.4,
		                                                            new string[] {"AAT", "AAC"});

		public static readonly AminoAcid AsparticAcid = new AminoAcid("C4H5NO3", "Aspartic Acid", "Asp", 'D', 5.9,
		                                                              new string[] {"GAT", "GAC"});

		public static readonly AminoAcid Cysteine = new AminoAcid("C3H5NOS", "Cysteine", "Cys", 'C', 3.3,
		                                                          new string[] {"TGT", "TGC"});

		public static readonly AminoAcid GlutamicAcid = new AminoAcid("C5H7NO3", "Glutamic Acid", "Glu", 'E', 5.8,
		                                                              new string[] {"GAA", "GAG"});

		public static readonly AminoAcid Glutamine = new AminoAcid("C5H8N2O2", "Glutamine", "Gln", 'Q', 3.7,
		                                                           new string[] {"CAA", "CAG"});

		public static readonly AminoAcid Glycine = new AminoAcid("C2H3NO", "Glycine", "Gly", 'G', 7.4,
		                                                         new string[] {"GGT", "GGC", "GGA", "GGG"});

		public static readonly AminoAcid Histidine = new AminoAcid("C6H7N3O", "Histidine", "His", 'H', 2.9,
		                                                           new string[] {"CAT", "CAC"});

		public static readonly AminoAcid Isoleucine = new AminoAcid("C6H11NO", "Isoleucine", "Ile", 'I', 3.8,
		                                                            new string[] {"ATT", "ATC", "ATA"});

		public static readonly AminoAcid Leucine = new AminoAcid("C6H11NO", "Leucine", "Leu", 'L', 7.6,
		                                                         new string[] {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"});

		public static readonly AminoAcid Lysine = new AminoAcid("C6H12N2O", "Lysine", "Lys", 'K', 7.2,
		                                                        new string[] {"AAA", "AAG"});

		public static readonly AminoAcid Methionine = new AminoAcid("C5H9NOS", "Methionine", "Met", 'M', 1.8,
		                                                            new string[] {"ATG"});

		public static readonly AminoAcid Phenylalanine = new AminoAcid("C9H9NO", "Phenylalanine", "Phe", 'F', 4,
		                                                               new string[] {"TTT", "TTC"});

		public static readonly AminoAcid Proline = new AminoAcid("C5H7NO", "Proline", "Pro", 'P', 5,
		                                                         new string[] {"CCT", "CCC", "CCA", "CCG"});

		public static readonly AminoAcid Serine = new AminoAcid("C3H5NO2", "Serine", "Ser", 'S', 8.1,
		                                                        new string[] {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"});

		public static readonly AminoAcid Threonine = new AminoAcid("C4H7NO2", "Threonine", "Thr", 'T', 6.2,
		                                                           new string[] {"ACT", "ACC", "ACA", "ACG"});

		public static readonly AminoAcid Tryptophan = new AminoAcid("C11H10N2O", "Tryptophan", "Trp", 'W', 1.3,
		                                                            new string[] {"TGG"});

		public static readonly AminoAcid Tyrosine = new AminoAcid("C9H9NO2", "Tyrosine", "Tyr", 'Y', 3.3,
		                                                          new string[] {"TAT", "TAC"});

		public static readonly AminoAcid Valine = new AminoAcid("C5H9NO", "Valine", "Val", 'V', 6.8,
		                                                        new string[] {"GTT", "GTC", "GTA", "GTG"});

		public static readonly AminoAcid[] aminoAcids =
			new AminoAcid[]
				{
					Alanine, Arginine, Asparagine, AsparticAcid,
					Cysteine, Glutamine, GlutamicAcid, Glycine,
					Histidine, Isoleucine, Leucine, Lysine,
					Methionine, Phenylalanine, Proline, Serine,
					Threonine, Tryptophan, Tyrosine, Valine
				};

		private static readonly Dictionary<SilacLabel, AminoAcid> labelList =
			new Dictionary<SilacLabel, AminoAcid>();

		public static readonly Dictionary<char, double> aaMonoMasses = InitMasses();

		public static readonly Dictionary<char, double> aaOccurences = InitOccurences();

		public static readonly Dictionary<char, double> aaWeights = InitWeights();

		public static readonly Dictionary<string, char> codonToAa = InitCodonToAa();

		public static readonly string singleLetterAas = ExtractSingleLetterAA();

		private readonly string abbreviation;
		private readonly string[] codons;
		private readonly List<Molecule> labelingDiffMolecule1 = new List<Molecule>();
		private readonly List<Molecule> labelingDiffMolecule2 = new List<Molecule>();
		private readonly List<string> labelingFormula = new List<string>();
		private readonly List<SilacLabel> labels = new List<SilacLabel>();
		private readonly char letter;

		private readonly double occurence;

		static AminoAcid() {
			string file = FileUtil.GetConfigPath() + "labels.txt";
			StreamReader sr = new StreamReader(file);
			string line = sr.ReadLine();
			string[] w = line.Split('\t');
			int compositionIndex = GetIndex(w, "Composition");
			if (compositionIndex == -1) {
				throw new Exception("No column named 'Composition' found in labels.txt.");
			}
			int mascotNameIndex = GetIndex(w, "Mascot Name");
			if (mascotNameIndex == -1) {
				throw new Exception("No column named 'Mascot Name' found in labels.txt.");
			}
			int shortNameIndex = GetIndex(w, "Short Name");
			if (shortNameIndex == -1) {
				throw new Exception("No column named 'Short Name' found in labels.txt.");
			}
			int aaIndex = GetIndex(w, "AA");
			if (aaIndex == -1) {
				throw new Exception("No column named 'AA' found in labels.txt.");
			}
			while ((line = sr.ReadLine()) != null) {
				if (line.Length == 0) {
					continue;
				}
				string[] ww = line.Split('\t');
				AddLabeling(ww[compositionIndex], ww[mascotNameIndex], ww[shortNameIndex], ww[aaIndex][0]);
			}
		}

		internal AminoAcid(string empiricalFormula, string name, string abbreviation, char letter, double occurence,
		                   string[] codons)
			: base(empiricalFormula, name) {
			this.abbreviation = abbreviation;
			this.letter = letter;
			this.occurence = occurence / 100.0;
			this.codons = codons;
		}

		public string Abbreviation {
			get { return abbreviation; }
		}

		public char Letter {
			get { return letter; }
		}

		//     IUB/GCG      Meaning     Complement   Staden/Sanger

		// A             A             T             A
		// C             C             G             C
		// G             G             C             G
		//T/U            T             A             T
		// M           A or C          K             M
		// R           A or G          Y             R
		// W           A or T          W             W
		// S           C or G          S             S
		// Y           C or T          R             Y
		// K           G or T          M             K
		// V        A or C or G        B             V
		// H        A or C or T        D             H
		// D        A or G or T        H             D
		// B        C or G or T        V             B
		//X/N     G or A or T or C    X/N            N
		//./~      gap character      ./~            -

		private static Dictionary<string, char> InitCodonToAa() {
			Dictionary<string, char> result = new Dictionary<string, char>();
			foreach (AminoAcid aa in aminoAcids) {
				foreach (string codon in aa.codons) {
					result.Add(codon, aa.letter);
				}
			}
			result.Add("TAG", '*');
			result.Add("TGA", '*');
			result.Add("TAA", '*');
			string[] keys = DataUtil.GetKeys(result);
			foreach (string key in keys) {
				if (key.Contains("T")) {
					result.Add(key.Replace('T', 'U'), result[key]);
				}
			}
			return result;
		}

		private static Dictionary<char, double> InitMasses() {
			Dictionary<char, double> result = new Dictionary<char, double>();
			foreach (AminoAcid aa in aminoAcids) {
				result.Add(aa.Letter, aa.GetMonoIsotopicMass());
			}
			result.Add('X', 111.000000);
			return result;
		}

		private static Dictionary<char, double> InitOccurences() {
			Dictionary<char, double> result = new Dictionary<char, double>();
			foreach (AminoAcid aa in aminoAcids) {
				result.Add(aa.Letter, aa.occurence);
			}
			result.Add('X', 1);
			return result;
		}

		private static Dictionary<char, double> InitWeights() {
			Dictionary<char, double> result = new Dictionary<char, double>();
			foreach (AminoAcid aa in aminoAcids) {
				result.Add(aa.Letter, aa.GetMolecularWeight());
			}
			return result;
		}

		private static int GetIndex(string[] w, string s) {
			for (int i = 0; i < w.Length; i++) {
				if (w[i].Equals(s)) {
					return i;
				}
			}
			return -1;
		}

		public static SilacLabel[] GetSilacLabels() {
			return DataUtil.GetKeys(labelList);
		}

		public static SilacLabel GetSilacLabelByShortName(string shortName) {
			foreach (SilacLabel l in DataUtil.GetKeys(labelList)) {
				if (l.ShortName.ToLower().Equals(shortName.ToLower())) {
					return l;
				}
			}
			return null;
		}

		private static void AddLabeling(string formula, string mascotName, string shortName, char aa) {
			int index = singleLetterAas.IndexOf(aa);
			if (index < 0) {
				throw new Exception("Illegal amino acid: " + aa);
			}
			aminoAcids[index].AddLabeling(formula, mascotName, shortName);
		}

		private void AddLabeling(string formula, string mascotName, string shortName) {
			labelingFormula.Add(formula);
			SilacLabel silacLabel = new SilacLabel(shortName, mascotName);
			labels.Add(silacLabel);
			Molecule[] x = GetDifferences(this, new Molecule(formula));
			labelingDiffMolecule1.Add(x[0]);
			labelingDiffMolecule2.Add(x[1]);
			labelList.Add(silacLabel, this);
		}

		public static string GetMascotModificationStringForLabel(SilacLabel label) {
			return label.MascotName;
		}

		private static string ExtractSingleLetterAA() {
			StringBuilder result = new StringBuilder();
			for (int i = 0; i < aminoAcids.Length; i++) {
				result.Append(aminoAcids[i].letter);
			}
			return result.ToString();
		}

		public static AminoAcid GetAminoAcidFromLabel(SilacLabel label) {
			if (labelList.ContainsKey(label)) {
				return labelList[label];
			}
			throw new Exception("Label is not in label list: " + label);
		}

		public double GetLabelingMassDiff(SilacLabel silacLabel) {
			return GetLabelingDiff2(silacLabel).GetMostLikelyMass(0.2)
			       - GetLabelingDiff1(silacLabel).GetMostLikelyMass(0.2);
		}

		public Molecule GetLabelingDiff1(SilacLabel silacLabel) {
			for (int i = 0; i < labels.Count; i++) {
				if (labels[i] == silacLabel) {
					return labelingDiffMolecule1[i];
				}
			}
			throw new Exception("Label is not in label list: " + silacLabel);
		}

		public Molecule GetLabelingDiff2(SilacLabel silacLabel) {
			for (int i = 0; i < labels.Count; i++) {
				if (labels[i] == silacLabel) {
					return labelingDiffMolecule2[i];
				}
			}
			throw new Exception("Label is not in label list: " + silacLabel);
		}

		public static AminoAcid[] FromLetters(char[] c) {
			AminoAcid[] result = new AminoAcid[c.Length];
			for (int i = 0; i < c.Length; i++) {
				result[i] = FromLetter(c[i]);
			}
			return result;
		}

		public static AminoAcid FromLetter(char c) {
			for (int i = 0; i < aminoAcids.Length; i++) {
				if (aminoAcids[i].Letter == c) {
					return aminoAcids[i];
				}
			}
			return null;
		}

		public SilacLabel GetFittingLabel(SilacLabel[] l) {
			foreach (SilacLabel sl in l) {
				if (GetAminoAcidFromLabel(sl).Equals(this)) {
					return sl;
				}
			}
			throw new Exception("No suitable label.");
		}

		public bool HasFittingLabel(SilacLabel[] l) {
			foreach (SilacLabel sl in l) {
				if (GetAminoAcidFromLabel(sl).Equals(this)) {
					return true;
				}
			}
			return false;
		}

		public static char[] GetSingleLetters(AminoAcid[] aas) {
			char[] result = new char[aas.Length];
			for (int i = 0; i < aas.Length; i++) {
				result[i] = aas[i].letter;
			}
			return result;
		}

		public override bool Equals(object obj) {
			if (this == obj) {
				return true;
			}
			if (obj is AminoAcid) {
				AminoAcid other = (AminoAcid) obj;
				return other.letter == letter;
			}
			return false;
		}

		public override int GetHashCode() {
			return letter + 1;
		}

		public static double CalcMolecularWeight(string sequence) {
			double result = MolUtil.WeightNormalCTerminus + MolUtil.WeightNormalNTerminus;
			foreach (char aa in sequence) {
				if (aaWeights.ContainsKey(aa)) {
					result += aaWeights[aa];
				}
			}
			return result;
		}

		public static double CalcMonoisotopicMass(string sequence) {
			double result = MolUtil.MassNormalCTerminus + MolUtil.MassNormalNTerminus;
			foreach (char aa in sequence) {
				if (aaMonoMasses.ContainsKey(aa)) {
					result += aaMonoMasses[aa];
				}
			}
			return result;
		}

		public static bool ValidSequence(string sequence) {
			foreach (char x in sequence) {
				if (singleLetterAas.IndexOf(x) < 0) {
					return false;
				}
			}
			return true;
		}
	}
}