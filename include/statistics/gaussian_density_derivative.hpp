#ifndef GAUSSIAN_DENSITY_DERIVATIVE_HPP
#define GAUSSIAN_DENSITY_DERUVATIVE_HPP

#include <vector>

#include <gromacs/utility/real.h>


/*!
 *
 */
class GaussianDensityDerivative
{
    public:

        std::vector<real> estimateApprox(
                const std::vector<real> &sample,
                const std::vector<real> &eval);
        std::vector<real> estimateDirect(
                const std::vector<real> &sample,
                const std::vector<real> &eval);

        // setter functions:
        void setBandWidth(real bw);
        void setDerivOrder(unsigned int r);
        void setErrorBound(real eps);

//    private:

        // internal constants:
        const real SQRT2PI_ = std::sqrt(2.0*M_PI);

        // internal variables:
        real bw_;
        unsigned int r_;
        real eps_;

        std::vector<real> coefA_;

        // 
        real estimDirectAt(
                const std::vector<real> &sample,
                real eval);
        real estimApproxAt(
                const std::vector<real> &sample,
                real eval);

        // calculation of coefficients:
        std::vector<real> setupCoefA();
        std::vector<real> setupCoefARef();
        void setupCoefB();
        real setupCoefQ();

        // internal utilities:
        unsigned int factorial(unsigned int n);
};

#endif

