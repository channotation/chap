#ifndef WEIGHTED_KERNEL_DENSITY_ESTIMATOR_HPP
#define WEIGHTED_KERNEL_DENSITY_ESTIMATOR_HPP

#include <vector>

#include "geometry/spline_curve_1D.hpp"
#include "statistics/kernel_density_estimator.hpp"


/*!
 * \brief Nadaraya-Watson kernel smoother.
 *
 * Performs essentially the same job as KernelDensityEstimator, but instead
 * of merely returning an estimate for the density of a point set that is 
 * normalised such that it integrates to one, it returns a continues function
 * that aims at smoothly interpolating between the given function values 
 * (interpreted as weights at the sample values).
 */
class WeightedKernelDensityEstimator : public KernelDensityEstimator
{
    public:

        // estimation interface:
        SplineCurve1D estimate(
                std::vector<real> &samples,
                std::vector<real> &weights);


    private:

        // evaluation of weighted density:
        std::vector<real> calculateWeightedDensity(
                std::vector<real> &samples,
                std::vector<real> &weights,
                std::vector<real> &evalPoints);
};

#endif

