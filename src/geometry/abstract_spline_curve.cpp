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

#include "geometry/abstract_spline_curve.hpp"


/*!
 * Getter method for spline curve degree.
 */
int
AbstractSplineCurve::degree() const
{
    return degree_;
}


/*!
 * Getter method for the number of control points in this spline curve.
 */
int
AbstractSplineCurve::nCtrlPoints() const
{
    return nCtrlPoints_;
}


/*!
 * Getter method for number of knots.
 */
int
AbstractSplineCurve::nKnots() const
{
    return nKnots_;
}


/*!
 * Getting method for the knot vector.
 */
std::vector<real>
AbstractSplineCurve::knotVector() const
{
    return knots_;
}


/*!
 * Getter method that returns that set of unique knots, i.e. the knot vector
 * without the repeated knots at start and end.
 */
std::vector<real>
AbstractSplineCurve::uniqueKnots() const
{
    // extract unique knots from knot vector:
    std::vector<real> uniqueKnots(
            knots_.begin() + degree_,
            knots_.end() - degree_);

    // return unique knots:
    return uniqueKnots;
}


/*!
 * Function to shift the spline parameter by a given offset. Internally,
 * this function simply subtracts the second element in the given RVec from 
 * each knot. An RVec is used instead of a simple real to provide a hook for
 * extending this method to shifts in angular coordinate in the case of 
 * SplineCurve3D.
 */
void
AbstractSplineCurve::shift(const gmx::RVec &shift)
{
    for(auto it = knots_.begin(); it != knots_.end(); it++)
    {
        *it -= shift[SS];
    }
}


/*!
 *  Given an evaluation point t, this functions finds the interval index idx 
 *  such that knot[idx] <= t < knot[idx + 1], i.e. the interval index required
 *  as starting point for the de Boor recursion. Note that the case where t is 
 *  identical to the final knot is treated as a special case and the interval
 *  is found such that knot[idx] <= t <= knot[idx + 1], i.e. with inclusive 
 *  upper boundary.
 */
int
AbstractSplineCurve::findInterval(const real &evalPoint)
{
    // initialise index:
    int idx = -1;

    // check if eval point is identical to last knot:
    if( evalPoint == knots_.back() )
    {
        // find iterator to interval bound:
        std::vector<real>::iterator it_lo = std::lower_bound(knots_.begin(), 
                                                             knots_.end(), 
                                                             evalPoint);
        // calculate interval index:
        idx = it_lo - knots_.begin() - 1;
    }
    else
    {
        // find iterator to interval point:
        std::vector<real>::iterator it_hi = std::upper_bound(knots_.begin(), 
                                                             knots_.end(), 
                                                             evalPoint);

        // calculate interval index:
        idx = it_hi - knots_.begin() - 1;
    }

    // return interval index:
    return idx;
}

