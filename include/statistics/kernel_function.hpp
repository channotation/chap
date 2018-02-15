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


#ifndef KERNEL_FUNCTION_HPP
#define KERNEL_FUNCTION_HPP

#include <memory>

#include "gromacs/utility/real.h"


/*!
 * \brief Abstract class specifying an interface for kernel functions.
 *
 * This class specified an interface for implementing kernel functions to be
 * used within KernelDensityEstimator. 
 */
class AbstractKernelFunction
{
    public:

        /*!
         * Operator for evaluating the kernel function at a given argument. Note 
         * that this only returns the argument dependent part of the kernel 
         * expression, any constant factor is left out and can be obtained via the
         * normalisingFactor() method.
         */
        virtual real operator()(real x) = 0;

        /*!
         * Getter method that returns the normalising factor of the kernel 
         * function. This is split off from the evaluation operator so that the 
         * associated computation does not have to be carried out each time the
         * kernel is evaluated (i.e. the constant factor can be pulled out of a 
         * sum.
         */
        virtual real normalisingFactor() = 0;
};


/*!
 * Shorthand notation for pointer to kernel function objects.
 */
typedef std::unique_ptr<AbstractKernelFunction> KernelFunctionPointer;


/*!
 * Enum for selection of kernel functions.
 */
enum eKernelFunction {eKernelFunctionGaussian};


/*!
 * \brief Factory class for creation of kernel functions.
 *
 * This class implements a static method for creating kernel functions of
 * a specific type.
 */
class KernelFunctionFactory
{
    public:

        static KernelFunctionPointer create(
                const eKernelFunction kernelFunction);
};


/*!
 * \brief Gaussian kernel function.
 *
 * The Gaussian kernel is defined as:
 *
 * \f[
 *      K(x) = \frac{1}{\sqrt{2\pi}} \exp\left( -\frac{1}{2} x^2 \right)
 * \f]
 */
class GaussianKernelFunction : public AbstractKernelFunction
{
    public:

        // evaluates non-constant part of kernel:
        virtual real operator()(real x);

        // returns constant prefactor:
        virtual real normalisingFactor();
};

#endif

