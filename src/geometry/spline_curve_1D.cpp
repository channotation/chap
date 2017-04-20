#include <iostream>

#include "geometry/spline_curve_1D.hpp"


/*
 * Constructor.
 */
SplineCurve1D::SplineCurve1D(int degree,
                             std::vector<real> knotVector,
                             std::vector<real> ctrlPoints)
{
    nCtrlPoints_ = ctrlPoints.size();
    nKnots_ = knotVector.size();
    degree_ = degree;


    // ensure minimal number of control points for given degree:
    if( nCtrlPoints_ < degree_ + 1 )
    {
        std::cerr<<"ERROR: Need at least d + 1 control points!"<<std::endl;
        std::cerr<<"d = "<<degree_
                 <<" and n = "<<nCtrlPoints_
                 <<std::endl;
        std::abort();
    }

    // ensure minimal number of knots:
    if( nKnots_ < nCtrlPoints_ + degree_ + 1 )
    {
        std::cerr<<"ERROR: Need at least n + d + 1 knots!"<<std::endl;
        std::cerr<<"d = "<<degree_
                 <<" and n = "<<nCtrlPoints_
                 <<" and k = "<<nKnots_
                 <<std::endl;
        std::abort();
    }

    // assign knot vector and control points:
    knotVector_ = knotVector;
    ctrlPoints_ = ctrlPoints;
}


/*
 * Destructor.
 */
SplineCurve1D::~SplineCurve1D()
{

}


/*
 * Public interface for spline evaluation. This function takes an evaluation 
 * point as an argument and returns the spline's value at this point. The 
 * actual evaluation is handled by different functions and the method argument
 * specifies which of these (a naive sum over all basis splines or de Boor's
 * recursive algorithm) should be used.
 *
 * The function can also evaluate the spline curve's derivative at a point and
 * the derivOrder argument is used to specify which order of the derivative 
 * should be evaluated. If deriOrder != 0 (the zeroth derivative is the 
 * function itself), then the method argument is ignored.
 */
real
SplineCurve1D::evaluate(real &evalPoint, 
                        unsigned int derivOrder, 
                        eSplineEvalMethod method)
{
    // introduce return variable:
    real value = 0.0;

    // check for extrapolation:
    if( evalPoint < knotVector_.front() )
    {
        return linearExtrap(evalPoint, knotVector_.front(), derivOrder);
    }
    if( evalPoint > knotVector_.back() )
    {
        return linearExtrap(evalPoint, knotVector_.back(), derivOrder);
    }

    // evaluation of spline function or its derivatives?
    if( derivOrder == 0 )
    {
        // which method should be used?
        if( method == eSplineEvalNaive )
        {
            value = evaluateNaive(evalPoint);        
        }
        else if( method == eSplineEvalDeBoor )
        {
            value = evaluateDeBoor(evalPoint);
        }
        else
        {
            std::cerr<<"ERROR: Requested spline evaluation method not available!"<<std::endl;
            std::abort();
        }
    }
    else
    {
        value = evaluateDeriv(evalPoint, derivOrder);
    }

    // return evaluation result:
    return value;
}


/*
 * Evaluation interface conveniently defined as operator.
 */
real
SplineCurve1D::operator()(real &evalPoint, 
                          unsigned int derivOrder, 
                          eSplineEvalMethod method)
{
    // actual compuatation is handled by evaluate method:
    return evaluate(evalPoint, derivOrder, method);
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
SplineCurve1D::evaluateNaive(real &evalPoint)
{
    // initialise return value as zero:
    real value = 0.0;

    // loop over all control points:
    for(int i = 0; i < nCtrlPoints_; i++)
    {
        value += ctrlPoints_[i] * B_(knotVector_, degree_, i, evalPoint);
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
SplineCurve1D::evaluateDeBoor(real &evalPoint)
{
    // find interval such that t_i <= evalPoint < t_i+1:
    int idx = findInterval(evalPoint);

    // return result of de Boor recursion:
    return deBoorRecursion(degree_, idx, evalPoint, ctrlPoints_);
}


/*
 * This method evaluates the spline curve's derivative of a given order by 
 * simply calculating the sum of the basis spline derivatives weighted by the
 * control points.
 */
real 
SplineCurve1D::evaluateDeriv(real &evalPoint, unsigned int derivOrder)
{
    // initialise derivative value as zero:
    real value = 0.0;

    // sum derivative value over all control points:
    for(unsigned int i = 0; i < ctrlPoints_.size(); i++)
    {
        value += ctrlPoints_[i] * D_(knotVector_, degree_, i, evalPoint, derivOrder);
    }

    // return derivative value:
    return value;
}


/*
 * This method provides linear extrapolation for evaluation points outside the
 * knot vector range. 
 */
real
SplineCurve1D::linearExtrap(real &evalPoint, real &boundary, unsigned int derivOrder)
{
    // evaluate curve or derivative?
    if( derivOrder == 0 )
    {
        // calculate offset and slope:
        real offset = evaluate(boundary, 0, eSplineEvalDeBoor);
        real slope = evaluate(boundary, 1, eSplineEvalDeBoor);

        // return linearly extrapolated value:
        return offset + slope*(evalPoint - boundary);
    }
    else if( derivOrder == 1 )
    {
        // simply return endpoint slope:
        return evaluate(boundary, 1, eSplineEvalDeBoor);
    }
    else
    {
        // for linear extrapolation, higher order derivative is always zero:
        return 0.0;
    }
}

