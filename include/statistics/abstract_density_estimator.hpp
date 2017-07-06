#ifndef ABSRTACT_DENSITY_ESTIMATOR_HPP
#define ABSRTACT_DENSITY_ESTIMATOR_HPP

#include "gromacs/utility/real.h"

#include "geometry/spline_curve_1D.hpp"


/*
 *
 */
class AbstractDensityEstimator
{
    public:

        // density estimation interface:
        virtual SplineCurve1D estimate(
                std::vector<real> &samples) = 0;


        // setter functions for parameters:
        void setBw(real bw);

    private:

};

#endif

