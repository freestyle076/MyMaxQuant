/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
namespace MaxQuant.Mol {
	public class SequenceDatabase {
		/// <summary>
		/// Regular expression which describes how to parse the fasta sequence header to 
		/// obtain the accession number.
		/// </summary>
		private readonly string accessionParseRule;

		/// <summary>
		/// Regular expression which describes how to parse the fasta sequence header to 
		/// obtain the protein description.
		/// </summary>
		private readonly string descriptionParseRule;

		/// <summary>
		/// The name of this database e.g. as known by Mascot.
		/// </summary>
		private readonly string name;

		public SequenceDatabase(string name, string accessionParseRule, string descriptionParseRule) {
			this.name = name;
			this.accessionParseRule = accessionParseRule;
			this.descriptionParseRule = descriptionParseRule;
		}

		/// <summary>
		/// The name of this database e.g. as known by Mascot.
		/// </summary>
		public string Name {
			get { return name; }
		}

		/// <summary>
		/// Regular expression which describes how to parse the fasta sequence header to 
		/// obtain the accession number.
		/// </summary>
		public string AccessionParseRule {
			get { return accessionParseRule; }
		}

		/// <summary>
		/// Regular expression which describes how to parse the fasta sequence header to 
		/// obtain the protein description.
		/// </summary>
		public string DescriptionParseRule {
			get { return descriptionParseRule; }
		}
	}
}