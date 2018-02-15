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


#ifndef LINEAR_SPLINE_INTERP_1D
#define LINEAR_SPLINE_INTERP_1D

#include <vector>

#include "geometry/spline_curve_1D.hpp"


/*!
 * \brief Functor for performing linear interpolation in one dimension.
 *
 * This class can be used to linearly interpolate a set of function values 
 * \f$ f_i = f(x_i) \f$ defined at a set of evaluation points \f$ x_i \f$. It 
 * will return a spline curve of degree one which goes through all input 
 * points:
 *
 * \f[
 *      s(x_i) = f_i
 * \f]
 *
 * Internally, this is done by setting up a B-spline as
 *
 * \f[
 *      s(x) = \sum_{i=1}^N f_i B_{i,1}(x)
 * \f]
 *
 * where \f$ N \f$ is the number of evaluation points.
 *
 * As a linear spline, the resulting curve will not be smooth, but it will be
 * tight, i.e. it is not prone to overshoot where the curvature of \f$ f(x) \f$
 * is large. This property is useful when probability density  functions, the 
 * value of which may never be negative. It is therefore used in classes 
 * derived from AbstractDensityEstimator.
 */
class LinearSplineInterp1D
{
    public:

        // interpolation interface:
        SplineCurve1D interpolate(
                const std::vector<real> &x,
                const std::vector<real> &f);
        SplineCurve1D operator()(
                const std::vector<real> &x,
                const std::vector<real> &f);

};

#endif

