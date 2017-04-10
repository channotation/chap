#include <algorithm>
#include <iostream>

#include "geometry/abstract_spline_curve.hpp"




/*
 *  Given an evaluation point t, this functions finds the interval index idx 
 *  such that knot[idx] <= t < knot[idx + 1], i.e. the interval index required
 *  as starting point for the de Boor recursion. Note that the case where t is 
 *  identical to the final knot is treated as a special case and the interval
 *  is found such that knot[idx] <= t <= knot[idx + 1], i.e. with inclusive 
 *  upper boundary.
 */
int
AbstractSplineCurve::findInterval(real &evalPoint)
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


/*
 * The function executes the de Boor Recursion (also known as Cox-de-Boor 
 * algorithm) for evaluating a spline curve at a given evaluation point. 
 */
real
AbstractSplineCurve::deBoorRecursion(int r, 
                                     int i, 
                                     real &evalPoint,
                                     const std::vector<real> &ctrlCoefs)
{
    // has bottom of recursion been reached?
    if( r == 0 )
    { 
        // return control point:
        return ctrlCoefs[i];
    }
    else
    {
        // calculate prefactor denominator:
        real alpha = knotVector_[i + degree_ + 1 - r] - knotVector_[i];
        
        // handle zero division case:
        if( alpha < std::numeric_limits<real>::epsilon() )
        {
            // define 0/0 = 0:
            alpha = 0.0;
        }
        else
        {
            // calculate de Boor coefficient:
            alpha = (evalPoint - knotVector_[i]) / alpha;
        }

        // call recursion:
        return (1.0 - alpha) * deBoorRecursion(r - 1, i - 1, evalPoint, ctrlCoefs) 
                     + alpha * deBoorRecursion(r - 1, i, evalPoint, ctrlCoefs);
    }
}


