#include <iostream>

#include "geometry/spline_curve_3D.hpp"


/*
 * Constructor.
 */
SplineCurve3D::SplineCurve3D(int degree,
                             std::vector<real> knotVector,
                             std::vector<gmx::RVec> ctrlPoints)
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
SplineCurve3D::~SplineCurve3D()
{

}


/*
 * Public interface for spline curve evaluation. This function Will return the 
 * (vector-valued) value of the spline curve (or its derivative) at a given
 * evaluation point.
 */
gmx::RVec
SplineCurve3D::evaluate(real &evalPoint,
                        unsigned int derivOrder,
                        eSplineEvalMethod method)
{
    // initialise return value:
    gmx::RVec value(0.0, 0.0, 0.0);

    // create individual vectors for each dimension:
    // TODO: is it more efficient to do this in constructor?
    std::vector<real> ctrlCoefsX;
    std::vector<real> ctrlCoefsY;
    std::vector<real> ctrlCoefsZ;
    for(unsigned int i = 0; i < nCtrlPoints_; i++)
    {
        ctrlCoefsX.push_back(ctrlPoints_.at(i)[0]);
        ctrlCoefsY.push_back(ctrlPoints_.at(i)[1]);
        ctrlCoefsZ.push_back(ctrlPoints_.at(i)[2]);
    }

    // evaluate spline function in each dimension:
    real valueX = evaluateSplineFun(evalPoint, ctrlCoefsX, derivOrder, method);
    real valueY = evaluateSplineFun(evalPoint, ctrlCoefsY, derivOrder, method);
    real valueZ = evaluateSplineFun(evalPoint, ctrlCoefsZ, derivOrder, method);

    // return result in vectorial form:
    return gmx::RVec(valueX, valueY, valueZ);
}


/*
 * Public evalaution interface conveniently defined as an operator.
 */
gmx::RVec
SplineCurve3D::operator()(real &evalPoint,
                          unsigned int derivOrder,
                          eSplineEvalMethod method)
{
    return evaluate(evalPoint, derivOrder, method);
}

