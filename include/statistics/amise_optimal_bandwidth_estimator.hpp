#ifndef AMISE_OPTIMAL_BANDWIDTH_ESTIMATOR_HPP
#define AMISE_OPTIMAL_BANDWIDTH_ESTIMATOR_HPP

#include <cmath>
#include <functional>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

#include "statistics/gaussian_density_derivative.hpp"


/*!
 *
 */
class AmiseOptimalBandwidthEstimator
{
    friend class AmiseOptimalBandwidthEstimatorTest;
    FRIEND_TEST(AmiseOptimalBandwidthEstimatorTest, 
                AmiseOptimalBandwidthEstimatorApproximateDerivativeTest);

    public:
 
        // public interface for bandwidth estimation:
        real estimate(
                std::vector<real> &samples);

    private:
       
        // 
        GaussianDensityDerivative gdd_;

        // constants:
        const real SQRTPI_ = std::sqrt(M_PI);
        const real SQRT2PI_ = std::sqrt(2.0 * M_PI);

        // density derivative functionals:
        inline real functionalPhi6(real sigma);
        inline real functionalPhi8(real sigma);        
        real functionalPhi(
                const std::vector<real> &sample,
                real bw,
                const int deriv);
        real functionalPhiFast(
                const std::vector<real> &samples,
                real bw,
                int deriv);

        // bandwidth to be used in derivative estimation:
        real gammaFactor_;
        inline real gammaFactor(const real phi4, const real phi6); 
        real gamma(real bw);
    
        // implicit equation for omptimal bandwidth:
        real optimalBandwidthEquation(
                const real bw,
                const std::vector<real> &samples);
};

#endif

