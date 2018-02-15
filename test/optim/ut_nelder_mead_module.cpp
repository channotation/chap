// CHAP - The Channel Annotation Package
// 
// Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
// Stephen J. Tucker
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


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

