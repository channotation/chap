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
    real curveLength = this -> length();

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
            real miLength = length(prevOldParam, mi);

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

    // find interpolating spline between points spaces equally in arc length:
    SplineCurve3D newSpl = Interp(arcParams, 
                                  arcPoints, 
                                  eSplineInterpBoundaryHermite);

    // update this spline:
    this -> nCtrlPoints_ = newSpl.nCtrlPoints_;
    this -> nKnots_ = newSpl.nKnots_;
    this -> knotVector_ = newSpl.knotVector_;
    this -> ctrlPoints_ = newSpl.ctrlPoints_;
}


/*
 * Calculates the length along the arc of the curve between the two parameter
 * values a and b.
 *
 * TODO: might need a faster algorithm here?
 */
real
SplineCurve3D::length(real &lo, real &hi)
{ 
    // find intervals of evaluation points:
    int idxLo = findInterval(lo);
    int idxHi = findInterval(hi);

//    real tLo = knotVector_[idxLo];
//    real tHi = knotVector_[idxHi];


    // initialise length as zero:
    real length = 0.0;

    // first segment:
    real tLim = std::min(knotVector_[idxLo + 1], hi);
    length += arcLengthBoole(lo, tLim);

    // loop over intermediate spline segments:
    if( idxHi - idxLo > 1 )
    {
        for(unsigned int i = idxLo + 1; i < idxHi; i++)
        {
            // add length of current segment:
            length += arcLengthBoole(knotVector_[i], knotVector_[i + 1]);
         }
    }

    // final segment:
    length += arcLengthBoole(knotVector_[idxHi], hi);
 
    // return Boole's rule estimate of length:
    return length;
}


/*
 * Convenience function to calculate arc length of curve between first and last
 * support point.
 */
real 
SplineCurve3D::length()
{
    return length(knotVector_.front(), knotVector_.back());
}


/*
 * Returns the tangent vector at the given evaluation point. Simply a wrapper
 * around the general evaluation function requesting the first derivative at 
 * the given point.
 */
gmx::RVec
SplineCurve3D::tangentVec(real &evalPoint)
{
    return evaluate(evalPoint, 1, eSplineEvalDeBoor); 
}


/*
 * Returns the speed of the curve at the given evaluation point, where speed
 * refers to the magnitude of the tangent vector.
 */
real
SplineCurve3D::speed(real &evalPoint)
{
    // calculate tangent vector:
    gmx::RVec tangentVec = this -> tangentVec(evalPoint);

    // return magnitude of tangent vector:
    return std::sqrt(tangentVec[0]*tangentVec[0] + 
                     tangentVec[1]*tangentVec[1] +
                     tangentVec[2]*tangentVec[2]);
}


/*
 * Uses Newton-Cotes quadrature of curve speed to determine the length of the 
 * arc between two given parameter values. The specific quadrature rule applied
 * is Boole's rule.
 */
real
SplineCurve3D::arcLengthBoole(real &lo, real &hi)
{
    // determine intermediate evaluation points:
    real h = (hi - lo)/4.0;
    real t2 = lo + 1.0*h;
    real t3 = lo + 2.0*h;
    real t4 = lo + 3.0*h;

    // evaluate speed at support points: 
    real s1 = speed(lo);
    real s2 = speed(t2);
    real s3 = speed(t3);
    real s4 = speed(t4);
    real s5 = speed(hi);

    // evaluate Boole's law:
    return 2.0*h/45.0*(7.0*s1 + 32.0*s2 + 12.0*s3 + 32.0*s4 + 7.0*s5);
}

