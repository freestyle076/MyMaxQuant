/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using MaxQuant.Spec;

namespace MaxQuant.Spec {
	public interface IIdentifiedPeptideList {
		int Count { get;}
		void Delete();
		void Move(string path1);
		void Write();
		void FillingFinished();
		bool IsUniqueGroup(int index);
		string GetSequence(int index);
		double GetPep(int index, PeptideScoring scoring, int nsilac);
	}
}