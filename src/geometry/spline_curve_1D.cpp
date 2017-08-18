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
    knots_ = knotVector;
    ctrlPoints_ = ctrlPoints;
}



/*
 *
 */
SplineCurve1D::SplineCurve1D()
{
    
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
    // use constant extrapolation here:
    if( evalPoint <= knotVector_.front() )
    {
        if( derivOrder > 0 )
        {
            return 0.0;
        }
        else
        {
            return ctrlPoints_.front();
        }
    }
    if(  evalPoint >= knotVector_.back() )
    {
        if( derivOrder > 0 )
        {
            return 0.0;
        }
        else
        {
            return ctrlPoints_.back();
        }
    }

    // one-dimensional case is just evaluation of spline function:
    return evaluateSplineFun(evalPoint, ctrlPoints_, derivOrder, method);
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
 *
 */
real
SplineCurve1D::evaluate(const real &eval, unsigned int deriv)
{
    // interpolation or extrapolation:
    if( eval < knots_.front() || eval > knots_.back() )
    {
        return evaluateExternal(eval, deriv);
    }
    else
    {
        return evaluateInternal(eval, deriv);
    }
}


/*
 *
 */
real
SplineCurve1D::evaluateInternal(const real &eval, unsigned int deriv)
{
    // container for basis functions or derivatives:
    SparseBasis basis;

    // derivative required?
    if( deriv == 0 )
    {
        // evaluate B-spline basis:
        basis = B_(eval, knots_, degree_);
    }
    else
    {
        // evaluate B-spline basis derivatives:
        basis = B_(eval, knots_, degree_, deriv); 
    }
    
    // return value of spline curve (derivative) at given evalaution point:
    return computeLinearCombination(basis);
 
}


/*
 *
 */
real
SplineCurve1D::evaluateExternal(const real &eval, unsigned int deriv)
{
    // which boundary is extrapolation based on?
    real boundary;
    if( eval < knots_.front() )
    {
        boundary = knots_.front();
    }
    else
    {
        boundary = knots_.back();
    }

    // derivative required?
    if( deriv == 0 )
    {
        // return value of curve at boundary:
        SparseBasis basis = B_(boundary, knots_, degree_);
        return computeLinearCombination(basis);
    }
    else
    {
        // for linear extrapolation, second and higher order deriv are zero:
        return 0.0;
    }
}


/*
 *
 */
real
SplineCurve1D::computeLinearCombination(const SparseBasis &basis)
{
    real value = 0.0; 
    for(auto b : basis)
    {
        value += b.second * ctrlPoints_[b.first];
    }

    return value;
}


/*!
 *  Getter function for access to the spline curves control points.
 */
std::vector<real>
SplineCurve1D::ctrlPoints() const
{
    return ctrlPoints_;
}

