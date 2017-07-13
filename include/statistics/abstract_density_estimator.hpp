#ifndef ABSRTACT_DENSITY_ESTIMATOR_HPP
#define ABSRTACT_DENSITY_ESTIMATOR_HPP

#include "gromacs/utility/real.h"

#include "geometry/spline_curve_1D.hpp"


/*!
 * \brief Helper class for specifying parameters is in the classes derived from
 * AbstractDensityEstimator.
 *
 * This class is used to simplify the interface of the various density 
 * estimation classes. It internally maintains variables for all parameters 
 * that may be used by any of these as well as corresponding flags indicating
 * whether the value of a specific parameter has been set.
 *
 * It is the responsibility of the classes derived from 
 * AbstractDensityEstimator to ensure that a parameter has been properly set
 * before using it. This class does not perform any sanity checks on the 
 * parameter values.
 */
class DensityEstimationParameters
{
    public:

        // constructor and destructor:
        DensityEstimationParameters();

        // setter methods:
        void setBinWidth(real binWidth);

        // getter methods:
        real binWidth() const;
        bool binWidthIsSet() const;

    private:

        real binWidth_;
        bool binWidthIsSet_;
    
};


/*!
 * \brief Abstract interface class for density estimation.
 *
 * This class specifies the interface that density estimation classes need to
 * implement. This enables substitution different methods (such as histograms
 * and kernel density estimation for one another without the need to rewrite
 * the code using these classes.
 */
class AbstractDensityEstimator
{
    public:

        // density estimation interface:
        virtual SplineCurve1D estimate(
                std::vector<real> &samples) = 0;

        // setter function for parameters:
        virtual void setParameters(
                const DensityEstimationParameters &params) = 0;

    private:


};

#endif

