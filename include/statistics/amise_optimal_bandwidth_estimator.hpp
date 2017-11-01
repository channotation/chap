#ifndef AMISE_OPTIMAL_BANDWIDTH_ESTIMATOR_HPP
#define AMISE_OPTIMAL_BANDWIDTH_ESTIMATOR_HPP

#include <cmath>
#include <functional>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"


/*!
 *
 */
class AmiseOptimalBandwidthEstimator
{
    public:
 
        // public interface for bandwidth estimation:
        real estimate(
                const std::vector<real> &samples);

    private:
        
        // constants:
        const real SQRTPI_ = std::sqrt(M_PI);
        const real SQRT2PI_ = std::sqrt(2.0 * M_PI);

        // density derivative functionals:
        inline real functionalPhi6(real sigma);
        inline real functionalPhi8(real sigma);        
        inline real functionalPhi(
                const std::vector<real> &samples,
                const real bw,
                const int deriv);
        inline real functionalPhiFast(
                std::vector<real> &samples,
                const real bw,
                const int deriv);

        // bandwidth to be used in derivative estimation:
        real gammaFactor_;
        inline real gammaFactor(const real phi4, const real phi6); 
        real gamma(real bw);
    
        // implicit equation for omptimal bandwidth:
        real optimalBandwidthEquation(
                const real bw,
                const std::vector<real> &samples);


        std::vector<real> intervalCentres(
                const real bw);

        int truncationNumber(
                const real bw,
                const real ir,
                const real cr,
                const real epsPrime,
                const unsigned int deriv);


        //
//        real cutoffRadius(const real );
};


#endif

