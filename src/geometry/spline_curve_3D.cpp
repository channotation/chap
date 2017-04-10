#include <iostream>
#include <cstdlib>
#include <limits>

#include "geometry/spline_curve_3D.hpp"


/*
 *
 */
SplineCurve3D::SplineCurve3D(int degree,
                             std::vector<real> knotVector,
                             std::vector<gmx::RVec> ctrlPoints)
    : degree_(degree)
    , nKnots_(knotVector.size())
    , nCtrlPoints_(ctrlPoints.size())
    , knotVector_(knotVector)
    , ctrlPoints_(ctrlPoints)
{
    // sanity checks:
    if( degree_ <= 0 )
    {
        std::cout<<"ERROR: Polynomial degree must be positive integer!"<<std::endl;
        abort();
    }
    if( nKnots_ != nCtrlPoints_ + degree_ + 1 )
    {
        std::cout<<"ERROR: Must have n + d + 1 knots, where n is number of control points and d is degree!"<<std::endl;
        std::cout<<"n = "<<nCtrlPoints_<<std::endl;
        std::cout<<"d = "<<degree_<<std::endl;
        abort();
    }
}


/*
 *
 */
SplineCurve3D::~SplineCurve3D()
{

}


/*
 * Public interface for evaluation of spline curve.
 */
gmx::RVec 
SplineCurve3D::evaluate(const real &evalPoint, 
                        eSplineEvalMethod method)
{
    // variable for the resturn value:
    gmx::RVec value;

    if( method == eSplineEvalNaive )
    {
        value = evalNaive(evalPoint);
    }
    else if( method == eSplineEvalDeBoor )
    {
        value = evalDeBoor(evalPoint);
    }

//    std::cout<<"value = "<<value[0]<<"  "<<value[1]<<"  "<<value[2]<<std::endl;

    // return result:
    return value;
}


/*
 * The naive method evaluates a spline simple through its definition as the 
 * sum over control-point weighted basis spline functions. This is not 
 * efficient and mainly implemented for testing and verification purposes.
 */
gmx::RVec
SplineCurve3D::evalNaive(const real &evalPoint)
{
   // std::cout<<"Naive"<<std::endl;
    // initialise return value as null vector:
    gmx::RVec value(0.0, 0.0, 0.0);

 //   std::cout<<"nKnots = "<<nKnots_<<std::endl;
 //   std::cout<<"nCtrlP = "<<nCtrlPoints_<<std::endl;

    // loop over all basis functions:
    for(int i = 0; i < nCtrlPoints_; i++)
    {
        // get value of basis spline at this evaluation point:
        real basisSplineVal = B_(knotVector_, degree_, i, evalPoint);
/*
        std::cout<<"i = "<<i<<"  "
                 <<"bspl = "<<basisSplineVal<<"  "
                 <<"ctrlP_x = "<<ctrlPoints_[i][0]<<"  "
                 <<"ctrlP_y = "<<ctrlPoints_[i][1]<<"  "
                 <<"ctrlP_z = "<<ctrlPoints_[i][2]<<std::endl;
*/
        // multiply onto control pointL
        value[0] += ctrlPoints_[i][0] * basisSplineVal;
        value[1] += ctrlPoints_[i][1] * basisSplineVal;
        value[2] += ctrlPoints_[i][2] * basisSplineVal;
    }

//    std::cout<<"value_xyz = "<<value[0]<<"  "<<value[1]<<"  "<<value[2]<<std::endl;

    // return evaluation result:
//    value = {6, 6, 6};
    return value;
}


/*
 *
 */
gmx::RVec
SplineCurve3D::evalDeBoor(const real &evalPoint)
{
    // check for extrapolation case:
    if( evalPoint < knotVector_.front() || evalPoint > knotVector_.back() )
    {
        std::cout<<"ERROR: Extrapolation is not implemented!"<<std::endl;
        abort();
    }

    // find interval:
    std::vector<real>::iterator it;
    it = std::lower_bound(knotVector_.begin(), knotVector_.end(), evalPoint);


    int idx = it - knotVector_.begin();
    std::cout<<"idx = "<<idx<<"  "
             <<"t[idx] = "<<knotVector_[idx]<<"  "
             <<"t[idx + 1] = "<<knotVector_[idx + 1]<<"  "
             <<"evalPoint = "<<evalPoint<<"  "
             <<std::endl;

    

    return deBoorRecursion(degree_, degree_, idx, evalPoint);
}


gmx::RVec
SplineCurve3D::deBoorRecursion(int degree, int r, int i, const real &t)
{
    gmx::RVec value;

    // check if bottom of recursion has been reached:
    if( r == 0 )
    {
        // return control point:
        value[0] = ctrlPoints_[i][0];
        value[1] = ctrlPoints_[i][1];
        value[2] = ctrlPoints_[i][2];
    }
    else
    {
        // calcykate prefactor denominator:
        real alpha = knotVector_[i + degree + 1 - r] - knotVector_[i];

        if( alpha <= std::numeric_limits<real>::epsilon() )
        {
//            std::cout<<"!!! zero division !!!"<<std::endl;   
            alpha = 0.0;
        }
        else
        {
            // calculate prefactor:
            alpha = (t - knotVector_[i]) / alpha;
        }

        // perform left and right recursion:
        gmx::RVec lRec = deBoorRecursion(degree, r - 1, i - 1, t);
        gmx::RVec rRec = deBoorRecursion(degree, r -1, i, t);

        // set return value:
        value[0] = (1 - alpha)*lRec[0] + alpha*rRec[0];
        value[1] = (1 - alpha)*lRec[1] + alpha*rRec[1];
        value[2] = (1 - alpha)*lRec[2] + alpha*rRec[2];
    }

    // return vectorial result:
    return value;
}






/*
 * De Boor's method is an efficient algorithm for evaluating a spline curve 
 * that makes use of the local propertz of splines to reduce the number of
 * necessary function evaluations.
 
gmx::RVec
SplineCurve3D::evalDeBoor(const real &evalPoint)
{
//    std::cout<<"deBoor"<<std::endl;
    // variable for return value:
    gmx::RVec value(0.0, 0.0, 0.0);

    // checj if evalPoint is in knot vector range:
    if( evalPoint < knotVector_.front() || evalPoint > knotVector_.back() )
    {
        std::cout<<"ERROR: Spline extrapolation not implemented!"<<std::endl;
        abort();
    }
    
    // find interval
    int interval = -1;
    for(int i = 0; i < nKnots_ - 1; i++)
    {
        if( evalPoint >= knotVector_[i] && evalPoint < knotVector_[i + 1] )
        {
            interval = i;
        }
        else
        {
            continue;
        }
    }
    if( interval == -1 )
    {
        interval = nKnots_ - 1;
    }

    

    for(int i = degree_ + 1; i <= nCtrlPoints_; i++)
    {
        // calculate zeroth order basis spline:
        real basisSplineVal = B_(knotVector_, 0, i, evalPoint);

        // calculate de Boor coefficients:
        gmx::RVec deBoorCoefs = deBoorRecursion(i, degree_, degree_, evalPoint);        
 
        // execute summation:
        value[0] += basisSplineVal * deBoorCoefs[0];
        value[1] += basisSplineVal * deBoorCoefs[1];
        value[2] += basisSplineVal * deBoorCoefs[2];
    }
  //  std::cout<<"value_xyz = "<<value[0]<<"  "<<value[1]<<"  "<<value[2]<<std::endl;

    // return evaluation result:
//    value = {9, 9, 9};
    return value;
}


gmx::RVec
SplineCurve3D::deBoorRecursion(int i, int d, int r,  const real &evalPoint)
{
    // check if bottom of recursion has been reached:
    if( d == 0 )
    {
        // simply return i-th control point:
        return ctrlPoints_[i];
    }
    else
    {
        // recursion index:
//        int r = d + 1 - i;
        std::cout<<"r = "<<r<<"  "
                 <<"d + 1 - i = "<<d + 1 - i<<std::endl;

        // calculate donominators:
        real denom = (knotVector_[i + r] - knotVector_[i]);
  
        real facA;
        real facB;

        // handle zero division case:
        if( denom <= std::numeric_limits<real>::epsilon() )
        {
            facA = 0.0;
            facB = 0.0;
        }
        else
        {
            facA = (knotVector_[i + r] - evalPoint) / denom;
            facB = (evalPoint - knotVector_[i]) / denom;
        }

        
        std::cout<<"r = "<<r<<"  "
                 <<"denom = "<<denom<<"  "
                 <<"facA = "<<facA<<"  "
                 <<"facB = "<<facB<<"  "
                 <<std::endl;
        
        // calculate prefactors:
        // calculate return value:
        gmx::RVec value;
        value[0] = facA*deBoorRecursion(i - 1, d, r - 1, evalPoint)[0] + facB*deBoorRecursion(i, d, r - 1, evalPoint)[0];
        value[1] = facA*deBoorRecursion(i - 1, d, r - 1, evalPoint)[1] + facB*deBoorRecursion(i, d, r - 1, evalPoint)[1];
        value[2] = facA*deBoorRecursion(i - 1, d, r - 1, evalPoint)[2] + facB*deBoorRecursion(i, d, r - 1, evalPoint)[2];

        // return value:
        return value;
    }
}
*/
