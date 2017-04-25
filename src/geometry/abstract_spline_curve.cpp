#include <algorithm>
#include <cmath>
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


/*
 *
 */
real
AbstractSplineCurve::evaluateSplineFun(real &evalPoint,
                                       const std::vector<real> &ctrlCoefs,
                                       unsigned int derivOrder, 
                                       eSplineEvalMethod method)
{
    // introduce return variable:
    real value = 0.0;

    // check for extrapolation:
    if( evalPoint < knotVector_.front() )
    {
        return linearExtrap(evalPoint, ctrlCoefs, knotVector_.front(), derivOrder);
    }
    if( evalPoint > knotVector_.back() )
    {
        return linearExtrap(evalPoint, ctrlCoefs, knotVector_.back(), derivOrder);
    }

    // evaluation of spline function or its derivatives?
    if( derivOrder == 0 )
    {
        // which method should be used?
        if( method == eSplineEvalNaive )
        {
            value = evaluateNaive(evalPoint, ctrlCoefs);        
        }
        else if( method == eSplineEvalDeBoor )
        {
            value = evaluateDeBoor(evalPoint, ctrlCoefs);
        }
        else
        {
            std::cerr<<"ERROR: Requested spline evaluation method not available!"<<std::endl;
            std::abort();
        }
    }
    else
    {
        value = evaluateDeriv(evalPoint, ctrlCoefs, derivOrder);
    }

    // return evaluation result:
    return value;
}


/*
 * Naive method for evaluating spline curve at given point that is based on
 * simply summing up all basis splines, i.e.
 *
 *     S(t) = sum_i^N C[i]*B(i, t)
 *
 * where C[i] is the i-th control point and B(t, i) the i-th basis spline 
 * evaluated at point t. 
 */
real
AbstractSplineCurve::evaluateNaive(real &evalPoint,
                                   const std::vector<real> &ctrlCoefs)
{
    // initialise return value as zero:
    real value = 0.0;

    // loop over all control points:
    for(int i = 0; i < nCtrlPoints_; i++)
    {
        value += ctrlCoefs[i] * B_(knotVector_, degree_, i, evalPoint);
    }

    // return evaluation result:
    return value;
}


/*
 * Efficient method for evaluating spline curve at given point that makes use
 * of de Boor's algorithm (and hence of the local support of the B-spline 
 * basis. This function is simply the driver for the de Boor recursion that 
 * finds the correct knot vector interval. The actual recursion is handled by
 * a separate function.
 */
real
AbstractSplineCurve::evaluateDeBoor(real &evalPoint,
                                    const std::vector<real> &ctrlCoefs)
{
    // find interval such that t_i <= evalPoint < t_i+1:
    int idx = findInterval(evalPoint);

    // return result of de Boor recursion:
    return deBoorRecursion(degree_, idx, evalPoint, ctrlCoefs);
}


/*
 * This method evaluates the spline curve's derivative of a given order by 
 * simply calculating the sum of the basis spline derivatives weighted by the
 * control points.
 */
real 
AbstractSplineCurve::evaluateDeriv(real &evalPoint, 
                                   const std::vector<real> &ctrlCoefs,
                                   unsigned int derivOrder)
{
    // initialise derivative value as zero:
    real value = 0.0;

    // find index to lowest interval where basis spline is nonzero:
    int idx = findInterval(evalPoint);

    // sum derivative value over all control points:
    for(int i = std::max(0, idx - degree_); i < std::min(nCtrlPoints_, idx + degree_); i++)
    {
        value += ctrlCoefs[i] * D_(knotVector_, degree_, i, evalPoint, derivOrder);
    }

    // return derivative value:
    return value;
}


/*
 * This method provides linear extrapolation for evaluation points outside the
 * knot vector range. 
 */
real
AbstractSplineCurve::linearExtrap(real &evalPoint, 
                                  const std::vector<real> ctrlCoefs,
                                  real &boundary, 
                                  unsigned int derivOrder)
{
    // evaluate curve or derivative?
    if( derivOrder == 0 )
    {
        // calculate offset and slope:
        real offset = evaluateSplineFun(boundary, ctrlCoefs, 0, eSplineEvalDeBoor);
        real slope = evaluateSplineFun(boundary, ctrlCoefs, 1, eSplineEvalDeBoor);

        // return linearly extrapolated value:
        return offset + slope*(evalPoint - boundary);
    }
    else if( derivOrder == 1 )
    {
        // simply return endpoint slope:
        return evaluateSplineFun(boundary, ctrlCoefs, 1, eSplineEvalDeBoor);
    }
    else
    {
        // for linear extrapolation, higher order derivative is always zero:
        return 0.0;
    }
}

