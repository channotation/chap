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


#include "geometry/linear_spline_interp_1D.hpp"

#include <stdexcept>


/*!
 * Public interpolation interface. Returns a one-dimensional spline curve that
 * linearly interpolates the given set of function values and evaluation 
 * points. Input vectors should be of same length and ordered in the same way.
 */
SplineCurve1D
LinearSplineInterp1D::interpolate(
        const std::vector<real> &x,
        const std::vector<real> &f)
{
    // sanity check:
    if( x.size() != f.size() )
    {
        throw std::logic_error("Input vectors for interpolation must be of "
        "the same size!");
    }

    // degree one for linear interpolation:
    int splineDegree = 1;

    // need to duplicate endpoints to get knot vector:
    std::vector<real> knotVector = x;
    knotVector.push_back(x.back());
    knotVector.insert(knotVector.begin(), x.front());

    // create spline object:
    SplineCurve1D splc(
        splineDegree,
        knotVector,
        f);

    // return the interpolating spline curve:
    return splc;
}


/*!
 * Interface to for linear interpolation as an operator. Simply forwards the
 * arguments to interpolate().
 */
SplineCurve1D
LinearSplineInterp1D::operator()(
        const std::vector<real> &x,
        const std::vector<real> &f)
{
    return interpolate(x, f);
}

