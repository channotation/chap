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
	std::cout<<"blub"<<std::endl;

	std::cout<<"stateDim"<<simAnMod.getStateDim()<<std::endl;

	// check that parameters have been set properly:
	ASSERT_EQ(stateDim, simAnMod.getStateDim());
}




int main(int argc, char **argv) {

	// initialise Google testing framework:
	testing::InitGoogleTest(&argc, argv);

	// run all tests:
	return RUN_ALL_TESTS();
}

