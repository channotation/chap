#ifndef ABSRTACT_DENSITY_ESTIMATOR_HPP
#define ABSRTACT_DENSITY_ESTIMATOR_HPP

#include "gromacs/utility/real.h"

#include "geometry/spline_curve_1D.hpp"


/*
 *
 */
class DensityEstimatorParameters
{
    public:

        void setBinWidth(real binWidth);

        real binWidth() const;

    private:

        real binWidth_;
    
};


/*
 *
 */
class AbstractDensityEstimator
{
    public:

        // density estimation interface:
        virtual SplineCurve1D estimate(
                std::vector<real> &samples) = 0;


        // setter function for parameters:
        virtual void setParameters(
                const DensityEstimatorParameters &params) = 0;

    private:


};



#endif

