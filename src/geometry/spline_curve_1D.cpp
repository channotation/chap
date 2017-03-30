#include <iostream>

#include "geometry/spline_curve_1D.hpp"


/*
 *
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
 *
 */
SplineCurve1D::~SplineCurve1D()
{

}


/*
 *
 */
real
SplineCurve1D::evaluate(real &evalPoint, eSplineEvalMethod method)
{
    // introduce return variable:
    real value = 0.0;

    // handle endpoint cases:


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
 *
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
 *
 */
real
SplineCurve1D::evaluateDeBoor(real &evalPoint)
{
    std::cout<<"Spline de Boor!"<<std::endl;

    // find interval such that t_i <= evalPoint < t_i+1:
    int idx = findInterval(evalPoint);

    // return result of de Boor recursion:
    return deBoorRecursion(degree_, idx, evalPoint, ctrlPoints_);
}

