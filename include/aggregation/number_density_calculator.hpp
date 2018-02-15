#ifndef NUMBER_DENSITY_CALCULATOR_HPP
#define NUMBER_DENSITY_CALCULATOR_HPP

#include <vector>

#include "gromacs/utility/real.h"

#include "geometry/spline_curve_1D.hpp"


/*!
 * \brief Functor for calculating physical number density.
 *
 * This class is used to convert the probability density \f$ p(s) \f$ of 
 * solvent particles along the pathway coordinate \f$ s \f$ into a physical 
 * number density by relating it to the local cross-sectional area according
 * to
 *
 * \f[
 *      n(s) = N \frac{p(s)}{\pi r(s)^2}
 * \f]
 *
 * where \f$ r(s) \f$ is the pathway radius at \f$ s \f$ and N is the overall
 * number of particles mapped onto the pathway. In other words, \f$ n(s) \f$
 * is the density of solvent particles in the reference volume
 *
 * \f[
 *      \pi r(s)^2 ds
 * \f]
 *
 * Note that as a one-dimensional probability density, \f$ p(s) \f$ has units 
 * of \f$ nm^{-1} \f$ so that \f$ n(s) \f$ has the appropriate units of 
 * \f$ nm^{-3} \f$. Multiplication with the overall number of particles yields
 * a physical number density that can be compared to number densities from
 * e.g. experiments. 
 */
class NumberDensityCalculator
{
    public:

        // interface for calculating number density:
        std::vector<real> operator()(
                const std::vector<real> &probabilityDensity,
                const std::vector<real> &radius,
                int totalNumber);

        
        // spline curve based version of interface:
        SplineCurve1D operator()(
                SplineCurve1D &probabilityDensity,
                SplineCurve1D &radius,
                int totalNumber);


    private:

        // auxiliary functions:
        std::vector<real> radiusToArea(
                const std::vector<real> &radius);
};

#endif

