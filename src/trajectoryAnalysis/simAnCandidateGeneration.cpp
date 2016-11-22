#include <cblas.h>

#include "trajectoryAnalysis/simAnCandidateGeneration.hpp"

/*
 *
 */
void
isotropicCandidateGeneration(int &dim,
                             real *crntState,
			     			 real *candState,
                             real *stateDir,
			     			 real *mat)
{
	// copy current state into candidate state:
	cblas_scopy(dim, crntState, 1, candState, 1);
	
	// add direction vector to candidate state:
	cblas_saxpy(dim, 1.0, stateDir, 1, candState, 1);
}


/*
 *
 */
void 
adaptiveCandidateGeneration(int &dim,
							real *crntState,
							real *candState,
							real *stateDir,
							real *mat)
{
	// copy current state into candidate state:
	cblas_scopy(dim, crntState, 1, candState, 1);

	// add adaptive direction to candidate state:
	cblas_sgemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, mat, 1, stateDir, 
				1, 1.0, candState, 1);
}
