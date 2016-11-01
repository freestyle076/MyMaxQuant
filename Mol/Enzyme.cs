/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Text;

namespace MaxQuant.Mol {
	public class Enzyme {
		private readonly CleavageRule[] cleavageRules;
		private readonly string name;

		public Enzyme(string name, string[] rules) {
			this.name = name;
			cleavageRules = new CleavageRule[rules.Length];
			for (int i = 0; i < rules.Length; i++) {
				cleavageRules[i] = new CleavageRule(rules[i]);
			}
		}

		public CleavageRule[] CleavageRules {
			get { return cleavageRules; }
		}

		public string Name {
			get { return name; }
		}

		public string CleavageSites {
			get {
				StringBuilder result = new StringBuilder();
				foreach (CleavageRule rule in cleavageRules) {
					result.Append(new string(rule.CleaveageSites.ToArray()));
				}
				return result.ToString();
			}
		}
	}
}