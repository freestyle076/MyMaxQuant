/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
namespace MaxQuant.Spec {
	public interface IIdentificationProvider {
		Identifications GetIdentifications(string rawfile, MascotQueryType type);
		void Dispose();
	}
}
