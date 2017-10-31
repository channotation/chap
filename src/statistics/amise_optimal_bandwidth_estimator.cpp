#include <boost/math/special_functions/factorials.hpp>

#include "statistics/amise_optimal_bandwidth_estimator.hpp"
#include "statistics/summary_statistics.hpp"


/*!
 *
 */
real
AmiseOptimalBandwidthEstimator::estimate(
        const std::vector<real> &samples)
{
    // sanity checks:
    if( samples.size() < 2 )
    {
        throw std::logic_error("Can not use AMISE optimal bandwidth "
                               "estimation with fewer than two samples.");
    }

    // 1. estimate standard deviation:
    SummaryStatistics sumStats;
    for(auto s : samples)
    {
        sumStats.update(s);
    }
    real sigma = sumStats.sd();

    // 2. estimate functionals six and eight from rule of thumb:
    real phi6 = functionalPhi6(sigma);
    real phi8 = functionalPhi8(sigma);
    

    // 3. estimate functionals four and six from KDE



    // 4. solve implicit bandwidth equation
}


/*!
 * Calculates the eighth order density derivative functional:
 *
 * \f[
 *      \phi_6 = \frac{-15}{16 \sqrt{\pi} \sigma^7}
 * \f]
 */
real
AmiseOptimalBandwidthEstimator::functionalPhi6(real sigma)
{
    return -15.0 / ( 16.0 * std::pow(sigma, 7) * SQRTPI_ );
}


/*!
 * Calculates the eighth order density derivative functional:
 *
 * \f[
 *      \phi_8 = \frac{105}{32 \sqrt{\pi} \sigma^9}
 * \f]
 */
real
AmiseOptimalBandwidthEstimator::functionalPhi8(real sigma)
{
    return 105.0 / ( 32.0 * std::pow(sigma, 9) * SQRTPI_ );
}

