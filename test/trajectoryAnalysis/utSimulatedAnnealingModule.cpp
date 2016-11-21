#include "trajectoryAnalysis/simulatedAnnealingModule.hpp"

#include <gtest/gtest.h>



TEST(utSimulatedAnnealingModule, constructorTest)
{
	// set parameters:
	int stateDim = 1;
	int seed = 15011991;
	int maxIter = 100;
	real temp = 310;
	real coolingFactor = 0.99;

	std::cout<<"blah"<<std::endl;

	// constructe a simulated annealing module:
	simulatedAnnealingModule simAnMod(stateDim, 
									  seed,
									  maxIter,
									  temp,
									  coolingFactor);

	// check that parameters have been set properly:
	ASSERT_EQ(stateDim, simAnMod.getStateDim());
	ASSERT_EQ(seed, simAnMod.getSeed());
	ASSERT_EQ(maxIter, simAnMod.getMaxIter());
	ASSERT_FLOAT_EQ(temp, simAnMod.getTemp());	// TODO: What is gromacs is compiled to use double?
	ASSERT_FLOAT_EQ(coolingFactor, simAnMod.getCoolingFactor());
}




int main(int argc, char **argv) {

	// initialise Google testing framework:
	testing::InitGoogleTest(&argc, argv);

	// run all tests:
	return RUN_ALL_TESTS();
}

