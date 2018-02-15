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


#include <iostream>

#include "geometry/spline_curve_1D.hpp"


/*!
 * Constructor for creating a spline curve of given degree from a set of knots
 * and control points.
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
    knots_ = knotVector;
    ctrlPoints_ = ctrlPoints;
}


/*!
 * Default constror for initialiser lists.
 */
SplineCurve1D::SplineCurve1D()
{
    
}


/*!
 * Public interface for evaluating the spline curve. In contrast to the method
 * of the same name in SplineCurve3D, this will use constant rather than linear
 * extrapolation when the evaluation point lies outside the range covered by
 * the knot vector.
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


/*!
 * Helper function for evaluating the spline curve at points inside the range 
 * covered by the knot vector.
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


/*!
 * Helper function for evaluating the spline curve at points outside the range
 * covered by the knot vector. Constant extrapolation is used here!
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


/*!
 * Auxiliary function for computing the linear combination of basis functions
 * weighted by control points.
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

