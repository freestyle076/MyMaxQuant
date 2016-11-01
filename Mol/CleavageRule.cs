/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using MaxQuant.Util;

namespace MaxQuant.Mol {
	public class CleavageRule {
		private readonly HashSet<char> cleaveageSites;
		private readonly bool cterm;
		private readonly HashSet<char> restrictions;

		public CleavageRule(string rule) {
			cterm = rule[0] == 'C';
			int tmp = rule.IndexOf("-");
			string cleaveSite;
			if (tmp == -1) {
				cleaveSite = rule.Substring(2);
				restrictions = new HashSet<char>();
			} else {
				cleaveSite = rule.Substring(2, tmp - 2);
				restrictions = new HashSet<char>();
				foreach (char cc in rule.Substring(tmp + 1)) {
					Restrictions.Add(cc);
				}
			}
			cleaveageSites = new HashSet<char>();
			foreach (char c in cleaveSite) {
				CleaveageSites.Add(c);
			}
		}

		public bool Cterm {
			get { return cterm; }
		}

		public HashSet<char> CleaveageSites {
			get { return cleaveageSites; }
		}

		public HashSet<char> Restrictions {
			get { return restrictions; }
		}
	}
}