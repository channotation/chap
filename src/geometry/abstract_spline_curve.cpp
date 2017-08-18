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
    return knotVector_;
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
            knotVector_.begin() + degree_,
            knotVector_.end() - degree_);

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
    for(auto it = knotVector_.begin(); it != knotVector_.end(); it++)
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
    if( evalPoint == knotVector_.back() )
    {
        // find iterator to interval bound:
        std::vector<real>::iterator it_lo = std::lower_bound(knotVector_.begin(), 
                                                             knotVector_.end(), 
                                                             evalPoint);
        // calculate interval index:
        idx = it_lo - knotVector_.begin() - 1;
    }
    else
    {
        // find iterator to interval point:
        std::vector<real>::iterator it_hi = std::upper_bound(knotVector_.begin(), 
                                                             knotVector_.end(), 
                                                             evalPoint);

        // calculate interval index:
        idx = it_hi - knotVector_.begin() - 1;
    }

    // return interval index:
    return idx;
}

