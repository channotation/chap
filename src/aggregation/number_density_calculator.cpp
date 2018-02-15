// CHAP - The Channel Annotation Package
// 
// Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
// Stephen J. Tucker
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "aggregation/number_density_calculator.hpp"


/*!
 * Public interface for calculating the number density of particles given 
 * probability density and local radius, as well as overall number of particles
 * in the sample. 
 *
 * Note that the input probability density is assumed to by normalised to one!
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

