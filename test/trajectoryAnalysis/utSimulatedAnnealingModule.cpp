#include <vector>
#include <limits>

#include "trajectoryAnalysis/simulatedAnnealingModule.hpp"

#include <gtest/gtest.h>


/*
 * Test fixture for units tests of simulated annealing class.
 */
class SimulatedAnnealingModuleTest : public ::testing::Test
{
	protected:

		const int randomSeed_ = 15011991;
		
		
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
 *
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
								 maxIter, 
								 initTemp, 
	                             coolingFactor, 
								 stepLengthFactor, 
								 initState,
								 rosenbrockFunction);
	
	// check if integer parameters have been set correctly:
	ASSERT_EQ(stateDim, sam.getStateDim());
	ASSERT_EQ(randomSeed_, sam.getSeed());
	ASSERT_EQ(maxIter, sam.getMaxIter());

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
 *
 */
TEST_F(SimulatedAnnealingModuleTest, CoolingTest)
{
	
}



/*
 *
 */
TEST_F(SimulatedAnnealingModuleTest, RosenbrockTest)
{
	// set tolerance for floating point comparison:
	real resTol = 1e-6;
	real errTol = 1e-3;

	// set parameters:
	int stateDim = 2;
	int maxIter = 5000;
	real temp = 300;
	real coolingFactor = 0.98;
	real stepLengthFactor = 0.01;

	// set initial state:
	real initState[stateDim] = {0.0, 0.0};

	// construct a simulated annealing module:
	SimulatedAnnealingModule sam(stateDim, 
								 randomSeed_,
								 maxIter,
								 temp,
								 coolingFactor,
								 stepLengthFactor,
								 initState,
								 rosenbrockFunction);

	// perform annealing:
	sam.anneal();

	// assert residual:
	real residual = sam.getBestCost() - 0.0;
	ASSERT_NEAR(residual, 0.0, resTol);

	// assert error:
	real error = 0.0;
	for(int i = 0; i < stateDim; i++)
	{
		error += (sam.getBestStateAt(i) - 1.0) * (sam.getBestStateAt(i) - 1.0);
	}
	error = std::sqrt(error);
	ASSERT_NEAR(error, 0.0, errTol);
}






/*
 * The Rosenbrock function is a standard test function for optimisation 
 * problems. Its minimum is at (a,a^2) where f(x,y) = 0. Note the non-canonical
 * negative sign in the return statement, as we are testing maximisation here.
 */
real rosenbrockFunction(std::vector<real> arg)
{
	// set internal parameters:
	real a = 1;
	real b = 100;

	// for legibility:
	real x = arg.at(0);
	real y = arg.at(1);

	// return value of Rosenbrock function at this point:
	return -(a - x)*(a - x) - b*(y - x*x)*(y - x*x);
}



/*
 * This test maximises the negative Rosenbrock function in two dimensions and 
 * asserts that the residual and error are below a certain tolerance threshold.
 */
TEST(utSimulatedAnnealingModule, rosenbrockTest)
{
/*
	// set tolerance for floating point comparison:
	real resTol = 1e-6;
	real errTol = 1e-3;

	// set parameters:
	int stateDim = 2;
	int seed = 15011991;
	int maxIter = 5000;
	real temp = 300;
	real coolingFactor = 0.98;
	real stepLengthFactor = 0.01;

	// set initial state:
	std::vector<real> initState = {0.0, 0.0};

	// construct a simulated annealing module:
	simulatedAnnealingModule simAnMod(stateDim, 
									  seed,
									  maxIter,
									  temp,
									  coolingFactor,
									  stepLengthFactor,
									  initState,
									  rosenbrockFunction);

	// perform annealing:
	simAnMod.anneal();

	// assert residual:
	real residual = simAnMod.getBestEnergy() - 0.0;
	ASSERT_NEAR(residual, 0.0, resTol);

	// assert error:
	std::vector<real> bestState = simAnMod.getBestState();
	real error = std::sqrt( (bestState.at(0) - 1.0)*(bestState.at(0) - 1.0) + (bestState.at(1) - 1.0)*(bestState.at(1) - 1.0) );
	ASSERT_NEAR(error, 0.0, errTol);
	*/
}

