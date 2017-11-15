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
    friend class AmiseOptimalBandwidthEstimatorTest;
    FRIEND_TEST(AmiseOptimalBandwidthEstimatorTest, 
                AmiseOptimalBandwidthEstimatorApproximateDerivativeTest);

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
        real functionalPhi(
                const std::vector<real> &samples,
                real bw,
                const int deriv);
 //       real functionalPhiFast(
   //             std::vector<real> &samples,
     //           const real bw,
       //         const int deriv);

        // bandwidth to be used in derivative estimation:
        real gammaFactor_;
        inline real gammaFactor(const real phi4, const real phi6); 
        real gamma(real bw);
    
        // implicit equation for omptimal bandwidth:
        real optimalBandwidthEquation(
                const real bw,
                const std::vector<real> &samples);


/*        std::vector<real> intervalCentres(
                const real bw);
        std::vector<size_t> nearestIntervalCentre(
                const std::vector<real> &centres,
                const std::vector<real> &samples);

        int truncationNumber(
                const real bw,
                const real ir,
                const real cr,
                const real epsPrime,
                const unsigned int deriv);
*/
       /* 
        real approximateDensityDerivative(
                const std::vector<real> &samples,
                const std::vector<real> &centres,
                const std::vector<size_t> &centreIdx,
                real eval,
                real bw,
                real rc,
                real ri,
                real epsPrime,
                real q,
                int deriv,
                unsigned int p);*/

  /*      
        real coefA(
                int s,
                int t,
                int deriv);
        real coefB(
                int k,
                int t,
                int l,
                const std::vector<real> &samples,
                const std::vector<real> &centres,
                const std::vector<size_t> &centreIdx,
                real bw,
                real coefQ);
        real coefC(
                int k,
                int s,
                int t,
                int deriv,
                real eval,
                real centre,
                real bw);
        real coefQ(
                real bw,
                int deriv,
                int numSamples);

*/
        //
  /*      real cutoffRadius(
                const real &ri,
                const real &bw,
                const real &derivFactorial,
                const real &epsPrime);*/
};


#endif

