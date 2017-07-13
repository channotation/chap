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

