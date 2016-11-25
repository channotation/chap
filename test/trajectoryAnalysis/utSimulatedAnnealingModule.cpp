#include <vector>
#include <limits>

#include <lapacke.h>

#include "trajectoryAnalysis/simulatedAnnealingModule.hpp"

#include <gtest/gtest.h>


/*
 * Test fixture for units tests of simulated annealing class.
 */
class SimulatedAnnealingModuleTest : public ::testing::Test
{
	protected:

		// default parameters:
		bool useAdaptiveCandidateGeneration_;

		int randomSeed_;
		int stateDim_;
		int maxCoolingIter_;
		int numCovSamples_;

		real initTemp_;
		real coolingFactor_;
		real stepLengthFactor_;
	

		// constructor:
		SimulatedAnnealingModuleTest()
		{
			useAdaptiveCandidateGeneration_ = false;

			randomSeed_ = 15011991;
			stateDim_ = 2;
			maxCoolingIter_ = 10000;
			numCovSamples_ = 50;
			
			initTemp_ = 300;
			coolingFactor_ = 0.99;
			stepLengthFactor_ = 0.01;
		}


		// functions to be used in tests:
		static real rosenbrockFunction(real *arg)
		{
			// set internal parameters:
			real a = 1;
			real b = 100;

			// for legibility:
			real x = arg[0];
			real y = arg[1];

			// return value of Rosenbrock function at this point:
			return -(a - x)*(a - x) - b*(y - x*x)*(y - x*x);
		}

};



/*
 * Tests the constructor of the simulated annealing class.
 */
TEST_F(SimulatedAnnealingModuleTest, ConstructorTest)
{
	// set parameters:
	int stateDim = 2;
	int maxIter = 1000;
	real initTemp = 300;
	real coolingFactor = 0.98;
	real stepLengthFactor = 0.01;
	real initState[stateDim] = {3.5, -2.0};

	// create simulated annealing module:
	SimulatedAnnealingModule sam(stateDim, 
								 randomSeed_, 
								 maxCoolingIter_,
								 numCovSamples_,
								 initTemp, 
	                             coolingFactor, 
								 stepLengthFactor, 
								 initState,
								 rosenbrockFunction,
								 useAdaptiveCandidateGeneration_);
	
	// check if integer parameters have been set correctly:
	ASSERT_EQ(stateDim, sam.getStateDim());
	ASSERT_EQ(randomSeed_, sam.getSeed());
	ASSERT_EQ(maxCoolingIter_, sam.getMaxCoolingIter());
	ASSERT_EQ(numCovSamples_, sam.getNumCovSamples());

	// check of real parameters have been set correctly:
	ASSERT_FLOAT_EQ(initTemp, sam.getTemp());
	ASSERT_FLOAT_EQ(coolingFactor, sam.getCoolingFactor());
	ASSERT_FLOAT_EQ(stepLengthFactor, sam.getStepLengthFactor());
	
	// check that state variables have been set correctly:
	for(int i = 0; i < stateDim; i++)
	{
		ASSERT_FLOAT_EQ(initState[i], sam.getCrntStateAt(i));
	}

	// check that correct energy is computed for initial states:
	ASSERT_FLOAT_EQ(rosenbrockFunction(initState), sam.getCrntCost());
	ASSERT_FLOAT_EQ(rosenbrockFunction(initState), sam.getCandCost());
	ASSERT_FLOAT_EQ(rosenbrockFunction(initState), sam.getBestCost());	
}



/*
 * Check that the cooling schedule works correctly.
 */
TEST_F(SimulatedAnnealingModuleTest, CoolingTest)
{
	// set parameters:
	int stateDim = 2;
	int maxIter = 1000;
	real initTemp = 300;
	real coolingFactor = 0.98;
	real stepLengthFactor = 0.01;
	real initState[stateDim] = {3.5, -2.0};

	// create simulated annealing module:
	SimulatedAnnealingModule sam(stateDim, 
								 randomSeed_, 
								 maxCoolingIter_,
								 numCovSamples_,
								 initTemp, 
	                             coolingFactor, 
								 stepLengthFactor, 
								 initState,
								 rosenbrockFunction,
								 useAdaptiveCandidateGeneration_);
	
	// perform cooling step:
	sam.cool();

	// check that temperature has decreased by correct amount:
	ASSERT_FLOAT_EQ(sam.getTemp(), initTemp*coolingFactor);

	// perform 9 more cooling steps:
	int additionalCoolingSteps = 0;
	for(int i = 0; i < additionalCoolingSteps; i++)
	{
		sam.cool();
	}

	// check that successive cooling also works:
	ASSERT_FLOAT_EQ(sam.getTemp(), initTemp*std::pow(coolingFactor, additionalCoolingSteps + 1));
}



/*
 * Tests that the simulated annealing module can successfully maximise the 
 * negative Rosenbrock function. This is done by asserting tolerance thresholds 
 * for the error and residual of the maximisation problem.
 */
TEST_F(SimulatedAnnealingModuleTest, RosenbrockTest)
{
	// set tolerance for floating point comparison:
	real resTol = 1e-6;
	real errTol = 1e-3;

	// set parameters:
	maxCoolingIter_ = 10000;
	numCovSamples_ = 1;
	real temp = 3000;
	real coolingFactor = 0.98;
	real stepLengthFactor = 0.01;

	// set initial state:
	real initState[stateDim_] = {0.0, 0.0};

	// construct a simulated annealing module:
	SimulatedAnnealingModule sam(stateDim_, 
								 randomSeed_,
								 maxCoolingIter_,
								 numCovSamples_,
								 temp,
								 coolingFactor,
								 stepLengthFactor,
								 initState,
								 rosenbrockFunction,
								 useAdaptiveCandidateGeneration_);

	// perform annealing:
	sam.anneal();

	// assert residual:
	real residual = sam.getBestCost() - 0.0;
	ASSERT_NEAR(residual, 0.0, resTol);

	// assert error:
	real error = 0.0;
	for(int i = 0; i < stateDim_; i++)
	{
		error += (sam.getBestStateAt(i) - 1.0) * (sam.getBestStateAt(i) - 1.0);
	}
	error = std::sqrt(error);
	ASSERT_NEAR(error, 0.0, errTol);
}



/*
 *
 */
TEST_F(SimulatedAnnealingModuleTest, test)
{
	real vec1[stateDim_] = {-1, 1};
	real vec2[stateDim_] = {-2, 2};
	real vec3[stateDim_] = {-3, 3};

	real mat[stateDim_ * 3] = {1, 2, 3,
	                           4, 5, 6};

	const int m = 2;
	const int n = 1;
	const int ldb = 3;
	const int lda = 1;


	LAPACKE_slacpy(LAPACK_ROW_MAJOR,
				   'A',
				   m,
				   n,
			   	   vec1,
			   	   lda,
			   	   &mat[0],
			   	   ldb);

	LAPACKE_slacpy(LAPACK_ROW_MAJOR,
				   'A',
				   m,
				   n,
			   	   vec2,
			   	   lda,
			   	   &mat[1],
			   	   ldb);

	LAPACKE_slacpy(LAPACK_ROW_MAJOR,
				   'A',
				   m,
				   n,
			   	   vec3,
			   	   lda,
			   	   &mat[2],
			   	   ldb);




	for(int i=0; i<2; i++)
	{
		for(int j=0; j<3; j++)
		{
			std::cout<<mat[j + 3*i]<<"  ";
		}
		std::cout<<std::endl;
	}


	std::cout<<"=============================================================="
	         <<std::endl
			 <<std::endl
			 <<std::endl;


	// set parameters:
	maxCoolingIter_ = 1;
	numCovSamples_ = 5;
	useAdaptiveCandidateGeneration_ = true;
	initTemp_ = 300.0;
	stepLengthFactor_ = 0.01;

	// set initial state:
	real initState[stateDim_] = {3.0, 2.0};

	// construct a simulated annealing module:
	SimulatedAnnealingModule sam(stateDim_, 
								 randomSeed_,
								 maxCoolingIter_,
								 numCovSamples_,
								 initTemp_,
								 coolingFactor_,
								 stepLengthFactor_,
								 initState,
								 rosenbrockFunction,
								 useAdaptiveCandidateGeneration_);

	// perform annealing:
	sam.anneal();


}



























