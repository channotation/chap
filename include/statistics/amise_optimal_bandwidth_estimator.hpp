#ifndef AMISE_OPTIMAL_BANDWIDTH_ESTIMATOR_HPP
#define AMISE_OPTIMAL_BANDWIDTH_ESTIMATOR_HPP

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"


/*!
 *
 */
class AmiseOptimalBandwidthEstimator
{
    friend class AmiseOptimalBandwidthEstimatorTest;
    FRIEND_TEST(
            AmiseOptimalBandwidthEstimatorTest, 
            AmiseOptimalBandwidthEstimatorGaussianDensityDerivativeFunctional);


    public:
    
        real estimate(
                const std::vector<real> &samples);

    private:
        
        // constants:
        const real SQRTPI_ = std::sqrt(M_PI);

        //
        inline real gaussianDensityDerivativeFunctional(
                const real sigma,
                const unsigned int deriv);

        // density derivative functionals:
        inline real functionalPhi6(real sigma);
        inline real functionalPhi8(real sigma);
};


#endif

