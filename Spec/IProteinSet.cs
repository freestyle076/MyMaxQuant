/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using MaxQuant.Mol;

namespace MaxQuant.Spec {
	public interface IProteinSet {
		string GetName(int i);
		Protein Get(int i);
		int GetLength(int i);
	}
}
