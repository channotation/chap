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

