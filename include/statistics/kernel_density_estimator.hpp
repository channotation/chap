#ifndef KERNEL_DENSITY_ESTIMATOR_HPP
#define KERNEL_DENSITY_ESTIMATOR_HPP

#include <vector>

#include "gromacs/utility/real.h"

#include "geometry/spline_curve_1D.hpp"

#include "statistics/abstract_density_estimator.hpp"


/*
 *
 */
class KernelDensityEstimator : public AbstractDensityEstimator
{
    public:
        
        // density estimation interface:
        virtual SplineCurve1D estimate(
                const std::vector<real> &samples);

    private:

};

#endif

