#ifndef WEIGHTED_KERNEL_DENSITY_ESTIMATOR_HPP
#define WEIGHTED_KERNEL_DENSITY_ESTIMATOR_HPP

#include <vector>

#include "geometry/spline_curve_1D.hpp"
#include "statistics/kernel_density_estimator.hpp"


/*
 *
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

