#include <gtest/gtest.h>

#include "optim/nelder_mead_module.hpp"


/*
 *
 */
class NelderMeadModuleTest : public ::testing::Test
{

    public:

        // Rosenbrock function as objective function:
        static real rosenbrock(std::vector<real> arg)
        {
            // internal parameters:
            real a = 1.0;
            real b = 100.0;

            // for legibility:
            real x = arg[0];
            real y = arg[1];

            // return value of Rosenbrock function:
            return -(a - x)*(a - x) - b*(y - x*x)*(y - x*x);
        };

        // sphere function as objective function:
        static real sphere(std::vector<real> arg)
        {
            // initialise result as zero:
            real res = 0.0;

            // sum squares of all components:
            for(int i = 0; i < arg.size(); i++)
            {
                res += arg[i]*arg[i];
            }

            // return result:
            return -res;
        };

};


/*
 *
 */
TEST_F(NelderMeadModuleTest, NelderMeadModuleSphereTest)
{
    // create Nelder-Mead optimisation module:
    NelderMeadModule nmm;

    // set objective function:
    nmm.setObjFun(sphere);

    // set optimisation parameters:
    std::map<std::string, real> params;
    params["maxIter"] = 100;
    params["initShift"] = 1.0;
    nmm.setParams(params);

    // prepare initial guess:
    std::vector<real> guess = {5.0, 6.0};
    nmm.setInitGuess(guess);

    // perform optimisation and extract result:
    nmm.optimise();
    OptimSpacePoint optim = nmm.getOptimPoint();

    // assert correct location of minimum:
    ASSERT_NEAR(0.0, optim.first[0], std::numeric_limits<real>::epsilon());
    ASSERT_NEAR(0.0, optim.first[1], std::numeric_limits<real>::epsilon());

    // assert correct value of minimum:
    ASSERT_NEAR(0.0, optim.second, std::numeric_limits<real>::epsilon());
}


/*
 *
 */
TEST_F(NelderMeadModuleTest, NelderMeadModuleRosenbrockTest)
{
    // create Nelder-Mead optimisation module:
    NelderMeadModule nmm;

    // set objective function:
    nmm.setObjFun(rosenbrock);

    // set optimisation parameters:
    std::map<std::string, real> params;
    params["maxIter"] = 200;
    params["initShift"] = 1.0;
    nmm.setParams(params);

    // prepare initial guess:
    std::vector<real> guess = {1.5, 1.5};
    nmm.setInitGuess(guess);

    // perform optimisation and extract result:
    nmm.optimise();
    OptimSpacePoint optim = nmm.getOptimPoint();

    // assert correct location of minimum:
    ASSERT_NEAR(1.0, optim.first[0], 10.0*std::numeric_limits<real>::epsilon());
    ASSERT_NEAR(1.0, optim.first[1], 10.0*std::numeric_limits<real>::epsilon());

    // assert correct value of minimum:
    ASSERT_NEAR(0.0, optim.second, std::numeric_limits<real>::epsilon());
}

