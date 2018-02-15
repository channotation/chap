#include <vector>
#include <limits>

#include "optim/simulated_annealing_module.hpp"

#include <gtest/gtest.h>


/*!
 * \brief Test fixture for units tests of simulated annealing class.
 */
class SimulatedAnnealingModuleTest : public ::testing::Test
{
		protected:

				// default parameters:
				bool useAdaptiveCandidateGeneration_;

				int randomSeed_;
				int stateDim_;
				int maxCoolingIter_;
				int numCostSamples_;

				real xi_;
				real convRelTol_;
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
						numCostSamples_ = 500;
						
						xi_ = 3.0f;
						convRelTol_ = 0.001;
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

        static real rosenbrock(std::vector<real> arg)
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


/*!
 * Tests that the simulated annealing module can successfully maximise the 
 * negative Rosenbrock function. This is done by asserting tolerance thresholds 
 * for the error and residual of the maximisation problem. This test uses the
 * isotropic, i.e. non-adaptive algorithm.
 */
TEST_F(SimulatedAnnealingModuleTest, IsotropicRosenbrockTest)
{
		// set tolerance for floating point comparison:
		real resTol = 1e-6;
		real errTol = 1e-3;

    // create a simulated annealing module:
    SimulatedAnnealingModule sam;

    // set parameters:
    std::map<std::string, real> params;
	params["saUseAdaptiveCandidateGeneration"] = 0;
    params["saRandomSeed"] = randomSeed_;
    params["saMaxCoolingIter"] = 100000;
    params["saNumCostSamples"] = 500;
    params["saInitTemp"] = 3000;
    params["saCoolingFactor"] = 0.99;
    params["saStepLengthFactor"] = 0.001;
    sam.setParams(params);
  
    // set initial state:
    std::vector<real> guess = {0.0, 0.0};
    sam.setInitGuess(guess);

    // set objective function:
    sam.setObjFun(rosenbrock);

    // perform optimisation:
    sam.optimise();

    // retrieve result:
    OptimSpacePoint res = sam.getOptimPoint();

    // assert correct cost:
    ASSERT_NEAR(0.0, res.second, resTol);

    // assert error:
    real error = 0.0;
    for(int i = 0; i < stateDim_; i++)
    {
            error += (res.first[i] - 1.0) * (res.first[i] - 1.0);
    }
    error = std::sqrt(error);
    ASSERT_NEAR(error, 0.0, errTol);

    // assert correct individual coordinates:
    ASSERT_NEAR(1.0, res.first[0], errTol);
    ASSERT_NEAR(1.0, res.first[1], errTol);
}

