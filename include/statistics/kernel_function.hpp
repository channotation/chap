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

