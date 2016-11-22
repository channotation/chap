#ifndef SIM_AN_CANDIDATE_GENERATION
#define SIM_AN_CANDIDATE_GENERATION

#include <gromacs/utility/real.h>

/*
 * Generates a candidate state by simply move in a random direction from the 
 * current state.
 */
void
isotropicCandidateGeneration(int &dim,
			     real *crntState,
			     real *candState,
                             real *stateDir,
			     real *mat);


/*
 * Generates candidate state by taking into account an estimate of the local 
 * shape of the objective function.
 */
//void 
//adaptiveCanidateGeneration();


#endif

