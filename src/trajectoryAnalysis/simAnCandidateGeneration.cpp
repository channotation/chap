#include "trajectoryAnalysis/simAnCandidateGeneration.hpp"

void
isotropicCandidateGeneration(int &dim,
                             real *crntState,
			     			 real *candState,
                             real *stateDir,
			     			 real *mat)
{
	// loop over elements:
	for(int i=0; i<dim; i++)
	{
		// perform simple vector addition:
		candState[i] = crntState[i] + stateDir[i];
	}
}


