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
 */
real
SplineCurve1D::evaluate(real &evalPoint, eSplineEvalMethod method)
{
    // introduce return variable:
    real value = 0.0;

    // check for extrapolation:
    if( evalPoint < knotVector_.front() || evalPoint > knotVector_.back() )
    {
        std::cerr<<"ERROR: Evaluation lies outside interpolation interval!"<<std::endl;
        std::cerr<<"Extrapolation is not yet implemented."<<std::endl;
        std::abort();
    }

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

    // return evaluation result:
    return value;
}


/*
 * Evaluation interface conveniently defined as operator.
 */
real
SplineCurve1D::operator()(real &evalPoint, eSplineEvalMethod method)
{
    // actual compuatation is handled by evaluate method:
    return evaluate(evalPoint, method);
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

