#include <gtest/gtest.h>

#include "optim/nelder_mead_module.hpp"


/*!
 * \brief Test fixture for the Nelder Mead optimisation module.
 *
 * Defines the objective functions used in the tests.
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
            for(size_t i = 0; i < arg.size(); i++)
            {
                res += arg[i]*arg[i];
            }

            // return result:
            return -res;
        };

};


/*!
 * Tests the Nelder Mead Module on the two-dimensional sphere function. Both
 * the error and the residual are checked and are required to converge to 
 * within the machine precision.
 */
TEST_F(NelderMeadModuleTest, NelderMeadModuleSphereTest)
{
    // create Nelder-Mead optimisation module:
    NelderMeadModule nmm;

    // set objective function:
    nmm.setObjFun(sphere);

    // set optimisation parameters:
    std::map<std::string, real> params;
    params["nmMaxIter"] = 100;
    params["nmInitShift"] = 1.0;
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


/*!
 * Tests the Nelder Mead Module on the two-dimensional Rosenbrock function. 
 * Both the error and the residual are checked and are required to converge to
 * within two times the machine precision.
 */
TEST_F(NelderMeadModuleTest, NelderMeadModuleRosenbrockTest)
{
    // create Nelder-Mead optimisation module:
    NelderMeadModule nmm;

    // set objective function:
    nmm.setObjFun(rosenbrock);

    // set optimisation parameters:
    std::map<std::string, real> params;
    params["nmMaxIter"] = 200;
    params["nmInitShift"] = 1.0;
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

