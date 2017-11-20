#ifndef AMISE_OPTIMAL_BANDWIDTH_ESTIMATOR_HPP
#define AMISE_OPTIMAL_BANDWIDTH_ESTIMATOR_HPP

#include <cmath>
#include <functional>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

#include "statistics/gaussian_density_derivative.hpp"


/*!
 * \brief Estimates the AMISE-optimal bandwidth for kernel density estimation.
 *
 * This class provides an estimate of the AMISE-optimal bandwidth
 *
 * \f[ 
 *      h = \frac{\Phi(K)^{\frac{1}{5}}}{m_2^{\frac{2}{5}} \Phi(p^{(r)})^{\frac{1}{5}} N^{\frac{1}{5}}}
 * \f]
 *
 * where \f$ m_2 \f$ is the second moment of the kernel function and
 *
 * \f[
 *      \Phi(f) = \int f(x)^2 dx
 * \f]
 *
 * is a funcational of the probability density (derivative). Its value is 
 * analytically known for the kernel itself, but needs to be estimated from the 
 * data via the procedure implemented in GaussianDensityDerivative for the
 * second derivative of the a priori unknown density of the data. 
 *
 * This implements the plug-in bandwidth selectior of Sheather and Jones 
 * (1991).
 *
 * \note The current implementation assumes a Gaussian kernel.
 */
class AmiseOptimalBandWidthEstimator
{
    friend class AmiseOptimalBandwidthEstimatorTest;
    FRIEND_TEST(AmiseOptimalBandwidthEstimatorTest, 
                AmiseOptimalBandwidthEstimatorApproximateDerivativeTest);

    public:
 
        // public interface for bandwidth estimation:
        real estimate(
                const std::vector<real> &samples);

    private:
       
        // 
        GaussianDensityDerivative gdd_;

        // constants:
        const real SQRTPI_ = std::sqrt(M_PI);
        const real SQRT2PI_ = std::sqrt(2.0 * M_PI);

        // density derivative functionals:
        inline real functionalPhi6(real sigma);
        inline real functionalPhi8(real sigma);        
        inline real functionalPhiFast(
                const std::vector<real> &samples,
                real bw,
                int deriv);

        // bandwidth to be used in derivative estimation:
        real gammaFactor_;
        inline real gammaFactor(const real phi4, const real phi6); 
        inline real gamma(real bw);
    
        // implicit equation for omptimal bandwidth:
        real optimalBandwidthEquation(
                const real bw,
                const std::vector<real> &samples);
};

#endif

