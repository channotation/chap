#include <algorithm>
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
    // number of uniformly spaced arc length parameters:
    int nNew = 10*nCtrlPoints_;

    // create lookup table for arc length at knots:
    prepareArcLengthTable();
   
    // determine uniform arc length spacing:
    real arcLenStep = arcLengthTable_.back() / (nNew - 1);

    // initialise new control points and parameters:
    std::vector<real> newParams;
    std::vector<gmx::RVec> newPoints;

    // loop over uniformly spaced arc length intervals:
    for(int i = 0; i < nNew; i++)
    {
        // calculate target arc length:
        real newParam = i*arcLenStep;
        newParams.push_back(newParam);

        // find parameter value corresponding to arc length value:
        real oldParam = arcLengthToParam(newParam); 

        //  evaluate spline to get new point:
        newPoints.push_back(this -> evaluate(oldParam, 0, eSplineEvalDeBoor));
    }
     
    // interpolate new points to get arc length parameterised curve:
    CubicSplineInterp3D Interp;
    SplineCurve3D newSpl = Interp(newParams, 
                                  newPoints, 
                                  eSplineInterpBoundaryHermite);

    // update own parameters:
    this -> knotVector_ = newSpl.knotVector_;
    this -> ctrlPoints_ = newSpl.ctrlPoints_;
    this -> nKnots_ = newSpl.nKnots_;
    this -> nCtrlPoints_ = newSpl.nCtrlPoints_;
    this -> arcLengthTableAvailable_ = false;
}


/*
 * Calculates the length along the arc of the curve between the two parameter
 * values a and b.
 */
real
SplineCurve3D::length(real &lo, real &hi)
{ 
    // initialise length as zero::
    real length = 0.0;   
  
    // find intervals of evaluation points:
    int idxLo = findInterval(lo);
    int idxHi = findInterval(hi);

    // add distance in endpoint intervals:
    if( idxHi == idxLo )
    {
        length += arcLengthBoole(lo, hi);
    }
    else
    {
        length += arcLengthBoole(lo, knotVector_[idxLo + 1]);
        length += arcLengthBoole(knotVector_[idxHi], hi);
    }

    // if necessary, loop over intermediate spline segments and sum up lengths:
    if( idxHi - idxLo > 1 )
    {
        for(unsigned int i = idxLo + 1; i < idxHi; i++)
        {
            // add length of current segment:
            length += arcLengthBoole(knotVector_[i], knotVector_[i + 1]);
         }
    }
   
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


/*
 * Prepares a lookup table that associates an arc length with each knot, where
 * the first knot is assigned an arc length of zero.
 */
void
SplineCurve3D::prepareArcLengthTable()
{
    // build arc length lookup table:
    for(unsigned int i = 0; i < knotVector_.size(); i++)
    {
        arcLengthTable_.push_back(length(knotVector_.front(), knotVector_.at(i)));
    }
 
    // set flag:
    arcLengthTableAvailable_ = true;
}


/*
 * Returns the parameter value (in the current parameterisation, typically 
 * chord length) that corresponds to a given value of arc length. This is done
 * via bisection.
 */
real
SplineCurve3D::arcLengthToParam(real &arcLength)
{
    int maxIter = 100;
    real absTol = 10*std::numeric_limits<real>::epsilon();

    // sanity check for arc length table:
    if( arcLengthTableAvailable_ != true )
    {
        prepareArcLengthTable();
    }

    // handle upper endpoint:
    if( arcLength == arcLengthTable_.back() )
    {
        return knotVector_.back();
    }

    // find appropriate interval:
    std::pair<std::vector<real>::iterator, std::vector<real>::iterator> bounds;
    bounds = std::equal_range(arcLengthTable_.begin(), arcLengthTable_.end(), arcLength);
    int idxLo = bounds.second - arcLengthTable_.begin() - 1;
    int idxHi = bounds.second - arcLengthTable_.begin();

    // initialise bisection interval and lower limit:
    real tLo = knotVector_[idxLo];
    real tHi = knotVector_[idxHi];
    real tMi = (tHi - tLo)/2.0;
    real tLi = tLo;
    
    // target arc length within this interval:
    real targetIntervalLength = arcLength - arcLengthTable_[idxLo];

    // bisection:
    int iter = 0;
    while( iter < maxIter )
    {
        // update interval midpoint:
        tMi = (tLo + tHi)/2.0;

        // evaluate length of arc:
        real intervalLength = arcLengthBoole(tLi, tMi);

        // calculate error:
        real signedError = intervalLength - targetIntervalLength;
        real error = std::abs(signedError);

        // break out of loop once error tolerance has been achieved:
        if( error < absTol )
        {
            break;
        }
        
        // determine new interval:
        if( signedError > 0.0 )
        {
            tHi = tMi;
        }
        else
        {
            tLo = tMi;
        }

        // increment loop counter:
        iter++;
    }

    // return bisection result:
    return tMi;
}

