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

#include "statistics/kernel_function.hpp"


/*!
 * Factory method for objects derived from AbstractKernelFunction. Takes the 
 * type of kernel function to be created as an argument and returns a unique
 * pointer to a new object of this type. Throws an exception if the requested
 * function is not implemented.
 */
KernelFunctionPointer
KernelFunctionFactory::create(const eKernelFunction kernelFunction)
{
    // create a unique pointer to a kernel function:
    KernelFunctionPointer kfp;

    // which kernel function should be used:
    if( kernelFunction == eKernelFunctionGaussian )
    {
        // create new Gaussian kernel function:
        kfp.reset(new GaussianKernelFunction());
    }
    else
    {
        throw std::runtime_error("Requested kernel function not available.");
    }

    // return pointer to kernel function:
    return kfp;
}


/*!
 * Evaluates the non-constant part of the Gaussian kernel:
 *
 * \f[
 *      \exp\left( -\frac{1}{2} x^2 \right)
 * \f]
 */
real
GaussianKernelFunction::operator()(real x)
{
    return std::exp( -0.5*x*x );
}


/*!
 * Returns the constant prefactor of a Gaussian kernel:
 *
 * \f[
 *      \frac{1}{\sqrt{2\pi}}
 * \f]
 */
real
GaussianKernelFunction::normalisingFactor()
{
    return 1.0/std::sqrt( 2.0*M_PI );
}

