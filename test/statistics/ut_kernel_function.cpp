#include <cmath>

#include <gtest/gtest.h>

#include "statistics/kernel_function.hpp"


/*!
 * Test fixture for kernel functions derived from AbstractKernelFunction.
 */
class KernelFunctionTest : public ::testing::Test
{

};


/*!
 * Test for the Gaussian kernel function. Asserts that the kernel fulfills the
 * defining properties of a kernel, namely non-negativity, symmetry, and 
 * normalisation.
 */
TEST_F(KernelFunctionTest, KernelFunctionGaussianTest)
{
    // tolerance for floating point comparison:
    real eps = std::numeric_limits<real>::epsilon();

    // create a Gaussian kernel function:
    KernelFunctionPointer Kernel = KernelFunctionFactory::create(
            eKernelFunctionGaussian);

    // create a sample of evaluation points:
    size_t numEvalPoints = 10000;
    real evalLo = -10.0;
    real evalHi = 10.0;
    real evalStep = (evalHi - evalLo) / (numEvalPoints - 1);
    std::vector<real> evalPoints;
    for(size_t i = 0; i < numEvalPoints; i++)
    {
        evalPoints.push_back(evalLo + i*evalStep);
    }

    // assert nonnegativity of kernel:
    for(auto eval : evalPoints)
    {
        ASSERT_LE(0.0, Kernel -> operator()(eval));
    }

    // assert symmetry of kernel:
    for(auto eval : evalPoints)
    {
        ASSERT_NEAR(
                Kernel -> operator()(eval),
                Kernel -> operator()(-eval),
                eps);
    }

    // assert normalisation of kernel:
    real integral = 0.0;
    for(eval : evalPoints)
    {
        integral += Kernel -> operator()(eval);
    }
    integral *= Kernel -> normalisingFactor();
    integral *= evalStep;
    ASSERT_NEAR(1.0, integral, std::sqrt(eps));
}

