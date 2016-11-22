#include "trajectoryAnalysis/simAnCandidateGeneration.hpp"

#include <gtest/gtest.h>


/*
 * This test checks if the isotropic candidate generation works without
 * problems. This is done by (i) mapping a state vector onto itself (i.e. zero
 * direction vector) and (ii) by mapping a state vector onto the null vector
 * (i.e. by using the negative of the state vector as direction vector).
 */
TEST(utSimAnCandidateGeneration, isotropicCandidateGenerationTest)
{
	// define state arrays:
	int dim = 3;
	real crntState[3] = {1.0, 2.0, 3.0};
	real candState[3];
	real stateDir1[3] = {0.0, 0.0, 0.0};
	real stateDir2[3] = {-crntState[0], -crntState[1], -crntState[2]};
	real *mat;		 	// not needed in isotropic state generation

	// zero difference move:
	isotropicCandidateGeneration(dim,
	                             crntState,
				     candState,
				     stateDir1,
				     mat);

	// assert that candidate and current state are similar:
	for(int i=0; i<dim; i++)
	{
		ASSERT_FLOAT_EQ(candState[i], crntState[i]);
	}

	// move state back to null vector:
	isotropicCandidateGeneration(dim,
	                             crntState,
				     candState,
				     stateDir2,
				     mat);
	
	// assert that candidate state is now null vector:
	for(int i=0; i<dim; i++)
	{
		ASSERT_FLOAT_EQ(candState[i], 0.0);
	}
}
