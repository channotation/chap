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


#include <cmath>

#include <gtest/gtest.h>

#include "statistics/kernel_function.hpp"


/*!
 * \brief Test fixture for kernel functions derived from AbstractKernelFunction.
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

    // assert non-negativity of kernel:
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
    for(auto eval : evalPoints)
    {
        integral += Kernel -> operator()(eval);
    }
    integral *= Kernel -> normalisingFactor();
    integral *= evalStep;
    ASSERT_NEAR(1.0, integral, std::sqrt(eps));
}

