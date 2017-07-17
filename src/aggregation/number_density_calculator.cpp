#include <algorithm>
#include <stdexcept>

#include "aggregation/number_density_calculator.hpp"


/*!
 * Public interface for calculating the number density of particles given 
 * probability density and local radius, as well as overall number of particles
 * in the sample.
 */
std::vector<real>
NumberDensityCalculator::operator()(
        const std::vector<real> &probabilityDensity,
        const std::vector<real> &radius,
        int totalNumber)
{
    // sanity checks:
    if( probabilityDensity.size() != radius.size() )
    {
        throw std::logic_error("Probability density and area input vectors "
                               "must be of equal size!");
    }

    // calculate cross-sectional area:
    std::vector<real> area = radiusToArea(radius);

    // calculate number density:
    std::vector<real> numberDensity = probabilityDensity;
    for(size_t i = 0; i < numberDensity.size(); i++)
    {
        numberDensity[i] *= totalNumber / area[i];
    }

    // return number density:
    return numberDensity;
}


/*!
 * Auxiliary function for calculating the circular cross-sectional area 
 * corresponding to a given radius:
 *
 * \f[
 *      A = \pi r^2
 * \f]
 */
std::vector<real>
NumberDensityCalculator::radiusToArea(
        const std::vector<real> &radius)
{
    // copy radiu into new vector:
    std::vector<real> area = radius;

    // calculate circular cross-sectional area using lambda expression:
    std::transform(
            area.begin(), 
            area.end(), 
            area.begin(), 
            [](real x) -> real{return x*x*M_PI;});

    // return area:
    return area;
}

