#include <vector>
#include <limits>

#include "trajectoryAnalysis/simulatedAnnealingModule.hpp"

#include <gtest/gtest.h>


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
 * This test makes checks if internal variables are initialised properly in the
 * constructor.
 */
TEST(utSimulatedAnnealingModule, constructorTest)
{
	// set parameters:
	int stateDim = 2;
	int seed = 15011991;
	int maxIter = 100;
	real temp = 310;
	real coolingFactor = 0.99;
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

	// check that parameters have been set properly:
	ASSERT_EQ(stateDim, simAnMod.getStateDim());
	ASSERT_EQ(seed, simAnMod.getSeed());
	ASSERT_EQ(maxIter, simAnMod.getMaxIter());
	ASSERT_FLOAT_EQ(temp, simAnMod.getTemp());	// TODO: What is gromacs is compiled to use double?
	ASSERT_FLOAT_EQ(coolingFactor, simAnMod.getCoolingFactor());

	// are state vectors of correct size:
	ASSERT_EQ(stateDim, simAnMod.getCrntState().size());
	ASSERT_EQ(stateDim, simAnMod.getCandState().size());
	ASSERT_EQ(stateDim, simAnMod.getBestState().size());

	// are state vectors properly initialised:
	for(int i = 0; i < stateDim; i++)
	{
		ASSERT_FLOAT_EQ(initState.at(i), simAnMod.getCrntState().at(i));
		ASSERT_FLOAT_EQ(initState.at(i), simAnMod.getCandState().at(i));
		ASSERT_FLOAT_EQ(initState.at(i), simAnMod.getBestState().at(i));
	}

	// is correct energy computed for initial states:
	ASSERT_FLOAT_EQ(rosenbrockFunction(initState), simAnMod.getCrntEnergy());	
}


/*
 *
 */
TEST(utSimulatedAnnealingModule, rosenbrockTest)
{
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
}



int main(int argc, char **argv) {

	// initialise Google testing framework:
	testing::InitGoogleTest(&argc, argv);

	// run all tests:
	return RUN_ALL_TESTS();
}

