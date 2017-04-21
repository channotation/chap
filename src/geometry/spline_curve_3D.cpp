#include <iostream>
#include <iomanip>
#include <limits>

#include "geometry/spline_curve_3D.hpp"
#include "geometry/cubic_spline_interp_3D.hpp"


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


/*
 * Change the internal representation of the curve such that it is 
 * parameterised in terms of arc length.
 */
void
SplineCurve3D::arcLengthParam()
{
    real absTol = 1e-2;
    int maxIter = 100;


    // calculate length of curve:
    real curveLength = this -> length(absTol);

    // number of uniformly spaced arc length parameters:
    int nSupport = nCtrlPoints_;

    // step length for uniform arc length parameterisation:
    real arcLengthStep = curveLength/(nSupport - 1);

    // vectors to hold data for lookup table:
    std::vector<real> arcParams;
    std::vector<gmx::RVec> arcPoints;
    std::vector<real> oldParams;

    // helper variables for temporary values:
    real prevArcParam = 0.0;
    real prevOldParam = knotVector_.front();

    // build lookup table:
    for(unsigned int i = 0; i < nSupport; i++)
    {
        // calculate uniformly spaced arc paramater points:
        real arcParam = i*arcLengthStep;
        
        // initialise bisection interval:
        real lo = prevOldParam;
        real hi = knotVector_.back();
        real mi = (lo + hi)/2.0;

        // find corresponding parameter value of old parameterisation:
        int iter = 0;
        while( iter < maxIter )
        {
            // calculate length on current interval:
            real miLength = length(prevOldParam, mi, 1e-3);

            // calculate error:
            real signedError = miLength - arcParam + prevArcParam;
            real error = std::abs(signedError);

            std::cout.precision(5); 
            std::cout<<"iter = "<<iter<<"  "
                     <<"arcParam = "<<arcParam<<"  "
                     <<"lo = "<<lo<<"  "
                     <<"hi = "<<hi<<"  "
                     <<"mi = "<<mi<<"  "
                     <<"miLength = "<<miLength<<"  "
                     <<"signedErr = "<<signedError<<"  "
                     <<"prevOldParam = "<<prevOldParam<<"  "
                     <<std::endl;
            

            // error tolerance reached?
            if( error < absTol )
            {
                break;
            }

            // update bisection interval:
            if( signedError < 0 )
            {
                lo = mi;
            }
            else
            {
                hi = mi;
            }
            mi = (lo + hi)/2.0;

            // increment loop counter:
            iter++;
        }

        // update 'previous iteration' values:
        prevArcParam = arcParam;
        prevOldParam = mi;

        // add values to lookup table:
        arcParams.push_back(arcParam);
        oldParams.push_back(mi);

        // evaluate spline curve at uniform arc length points:
        arcPoints.push_back(gmx::RVec(evaluate(mi, 0, eSplineEvalDeBoor)));
    }

    // interpolate between points uniformly spaced wrt arc length:
    CubicSplineInterp3D Interp;

    // TODO: add method to interpolation class that allows specifying 
    // parameterisation explictly!
   
}


/*
 * Calculates the length along the arc of the curve between the two parameter
 * values a and b.
 *
 * TODO: might need a faster algorithm here?
 */
real
SplineCurve3D::length(real &lo, real &hi, const real &absTol)
{   
    // initialse evaluation points as number of control points:
    int nEvalPoints = nCtrlPoints_;

    // initialise chord length, previous value a, and error:
    real chordLength;
    real prevChordLength = std::numeric_limits<real>::infinity();
    real error = 100.0;

    // iteratively calculate chord length until error is sufficiently small:
    int i = 0;
    while( error > absTol )
    {
        // reset chord length:
        chordLength = 0.0;

        // calculate parameter step:
        real paramStep = (hi - lo)/nEvalPoints; 

        // evaluate calculate chord length over evaluation points:
        gmx::RVec pointA = evaluate(lo, 0, eSplineEvalDeBoor);
        gmx::RVec pointB;
        for(unsigned int i = 1; i < nEvalPoints; i++)
        {
            // evaluate spline curve at current parameter value:
            real evalPoint = lo + i*paramStep;
            pointB = evaluate(evalPoint, 0, eSplineEvalDeBoor);

            // calculate Euclidean distance between points:
            real dX = pointB[0] - pointA[0];
            real dY = pointB[1] - pointA[1];
            real dZ = pointB[2] - pointA[2];
            real dist = std::sqrt(dX*dX + dY*dY + dZ*dZ);

            // add to chord length:
            chordLength += dist;

            // new point becomes old point:
            pointA = pointB;
        }

        // evaluate error:
        error = std::abs((chordLength - prevChordLength));
        prevChordLength = chordLength;

        // double the number of evaluation points:
        nEvalPoints *= 2;
    }    

    return chordLength;
}


/*
 * Convenience function to calculate arc length of curve between first and last
 * support point.
 */
real SplineCurve3D::length(const real &absTol)
{
    return length(knotVector_.front(), knotVector_.back(), absTol);
}


















