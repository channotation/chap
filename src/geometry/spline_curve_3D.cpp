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


#include <algorithm>
#include <limits>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>

#include "geometry/spline_curve_3D.hpp"
#include "geometry/cubic_spline_interp_3D.hpp"


/*!
 * Constructor for creating a spline curve of given degree from a knot vector
 * and associated control points.
 */
SplineCurve3D::SplineCurve3D(
        int degree,
        std::vector<real> knotVector,
        std::vector<gmx::RVec> ctrlPoints)
{
    // TODO: probably better to put some of these things into the initialiser
    // list

    nCtrlPoints_ = ctrlPoints.size();
    nKnots_ = knotVector.size();
    degree_ = degree;
    arcLengthTableAvailable_ = false;

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
 * Default constructor for initialiser lists. Does not set any members!
 */
SplineCurve3D::SplineCurve3D()
{

}


/*!
 * Public interface for evaluation of spline curve. Returns the value of the 
 * spline curve or its derivative at the given evaluation point. If the 
 * evaluation point lies outside the knot range, linear extrapolation is used.
 */
gmx::RVec
SplineCurve3D::evaluate(
        const real &eval,
        unsigned int deriv)
{
    // extrapolation or interpolation?
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
 * Auxiliary function for evaluating the spline curve at points inside the 
 * range covered by knots.
 */
gmx::RVec 
SplineCurve3D::evaluateInternal(const real &eval, unsigned int deriv)
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
    
    // return value of spline curve (derivative) at given evaluation point:
    return computeLinearCombination(basis);
}


/*!
 * Auxiliary function for evaluating the spline curve at points outside the 
 * range covered by knots. Linear extrapolation is used in this case.
 */
gmx::RVec 
SplineCurve3D::evaluateExternal(const real &eval, unsigned int deriv)
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
        // compute slope and offset:
        // TODO: this can be made more efficient by evaluating basis and derivs
        // in one go!
        SparseBasis basis = B_(boundary, knots_, degree_);
        gmx::RVec offset = computeLinearCombination(basis);
        basis = B_(boundary, knots_, degree_, 1);
        gmx::RVec slope = computeLinearCombination(basis);

        // return extrapolation point:
        svmul(eval - boundary, slope, slope);
        rvec_add(slope, offset, slope);
        return slope;

    }
    else if( deriv == 1 )
    {
        // simply return the slope at the endpoint:
        SparseBasis basis = B_(boundary, knots_, degree_, 1);
        return computeLinearCombination(basis);
    }
    else
    {
        // for linear extrapolation, second and higher order deriv are zero:
        return gmx::RVec(0.0, 0.0, 0.0);
    }
}


/*!
 * Evaluates the linear combination of basis functions weighted by control
 * points, i.e. computes
 *
 * \f[
 *      \sum_i \mathbf{c}_i B_{i,p}(x)
 * \f]
 *
 * where the sum only goes over the nonzero elements of the basis for 
 * efficiency.
 */
gmx::RVec
SplineCurve3D::computeLinearCombination(const SparseBasis &basis)
{
    gmx::RVec value(gmx::RVec(0.0, 0.0, 0.0)); 
    for(auto b : basis)
    {
        gmx::RVec tmp;
        svmul(b.second, ctrlPoints_[b.first], tmp);
        rvec_add(value, tmp, value);
    }

    return value;
}



/*!
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
        newPoints.push_back(this -> evaluate(oldParam, 0));
    }

    // interpolate new points to get arc length parameterised curve:
    CubicSplineInterp3D Interp;
    SplineCurve3D newSpl = Interp(newParams, 
                                  newPoints, 
                                  eSplineInterpBoundaryHermite);

    // update own parameters:
    this -> knots_ = newSpl.knots_;
    this -> ctrlPoints_ = newSpl.ctrlPoints_;
    this -> nKnots_ = newSpl.nKnots_;
    this -> nCtrlPoints_ = newSpl.nCtrlPoints_;
    this -> arcLengthTableAvailable_ = false;

    // reset reference points for mapping:
    refPoints_.clear();
}


/*!
 * Calculates the length along the arc of the curve between the two parameter
 * values a and b.
 */
real
SplineCurve3D::length(const real &lo, const real &hi)
{ 
    // do we need to form a lookup table:
    if( arcLengthTableAvailable_ == false || arcLengthTable_.size() == 0 )
    {
        prepareArcLengthTable(); 
    }

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
        length += arcLengthBoole(lo, knots_[idxLo + 1]);
        length += arcLengthBoole(knots_[idxHi], hi);
    }

    // if necessary, loop over intermediate spline segments and sum up lengths:
    if( idxHi - idxLo > 1 )
    {
        for(int i = idxLo + 1; i < idxHi; i++)
        {
            // add length of current segment from lookup table:
            length += arcLengthTable_[i + 1] - arcLengthTable_[i];
        }
    }

    return length;
}


/*!
 * Convenience function to calculate arc length of curve between first and last
 * support point.
 */
real 
SplineCurve3D::length()
{
    return length(knots_.front(), knots_.back());
}


/*!
 * Returns the tangent vector at the given evaluation point. Simply a wrapper
 * around the general evaluation function requesting the first derivative at 
 * the given point.
 */
gmx::RVec
SplineCurve3D::tangentVec(const real &eval)
{
    return evaluate(eval, 1); 
}

/*!
 * Returns the normal vector at the evaluation point.
 */
gmx::RVec
SplineCurve3D::normalVec(const real &eval)
{
    return evaluate(eval, 2);
}


/*!
 * Returns the speed of the curve at the given evaluation point, where speed
 * refers to the magnitude of the tangent vector.
 */
real
SplineCurve3D::speed(const real &eval)
{
    // return magnitude of tangent vector:
    return norm( this -> tangentVec(eval) );
}


/*!
 * Takes point in Cartesian coordinates and returns that points coordinates in
 * the curvilinear system defined by the spline curve. Return value is an RVec,
 * which contains the following information:
 *
 *      [0] - distance along the arc of the curve
 *      [1] - squared (!) distance from the curve at closest point
 *      [2] - angular coordinate; this is not yet implemented
 *
 * Note that this function assumes that the curve is parameterised by arc 
 * length!
 *
 * \todo implement angular coordinate!
 */
gmx::RVec 
SplineCurve3D::cartesianToCurvilinear(const gmx::RVec &cartPoint)
{
    // find index of interval containing closest point on spline curve:
    unsigned int idx = closestSplinePoint(cartPoint);

    // find closest point on this interval:
    gmx::RVec proj = projectionInInterval(
            cartPoint, 
            knots_[idx + degree_], 
            knots_[idx + degree_ + 1]);

    // check neighbouring knot intervals and extrapolate if necessary:
    gmx::RVec altProj;

    // next lower knot interval:
    if( idx == 0 )
    {
        altProj = projectionInExtrapRange(cartPoint, -1.0); 
    }
    else
    {
        altProj = projectionInInterval(
                cartPoint, 
                knots_[idx + degree_ - 1],
                knots_[idx + degree_]);
    }

    // does alternative projection give closer point:
    if( altProj[RR] < proj[RR] )
    {
        proj = altProj;
    }
    
    // next higher knot interval:
    if( idx == refPoints_.size() - 2 )
    {
        altProj = projectionInExtrapRange(cartPoint, 1.0); 
    }
    else
    {
        altProj = projectionInInterval(
                cartPoint, 
                knots_[idx + degree_ + 1],
                knots_[idx + degree_ + 2]);
    }

    // does alternative projection give closer point:
    if( altProj[RR] < proj[RR] )
    {
        proj = altProj;
    }
  
    // TODO: calculate angular coordinate!
    proj[PP] = 0.0;

    // return point in curvilinear coordinates:
    return proj;
}


/*!
 * Auxiliary function for finding the closest point on a spline curve that 
 * returns the corresponding spline interval index. First, a set of reference
 * points is sampled from the spline curve at the location of the unique knots
 * (this step is skipped if the reference points have already been computed in
 * a previous call to this function). Secondly, a linear search over these 
 * reference points is carried out to find the point closest to a given test 
 * point.
 *
 * The return value is the index of the closest reference point, except for the 
 * case where the closest reference point is the last point, which is mapped to 
 * the last interval, i.e. the index of the penultimate reference point is 
 * returned in this case.
 *
 * \todo Linear search may not be the most efficient.
 */
unsigned int
SplineCurve3D::closestSplinePoint(const gmx::RVec &point)
{
    // build lookup table:
    if( refPoints_.empty() )
    {
        refPoints_.reserve(uniqueKnots().size());
        for(auto s : uniqueKnots())
        {
            refPoints_.push_back( this -> evaluate(s, 0) );
        }
    }

    // find index of closest reference point on spline curve:
    unsigned int idxMinDist = 0;
    real minDist = std::numeric_limits<real>::infinity();
    for(unsigned int i = 0; i < refPoints_.size(); i++)
    {
        // find dist to control point:
        real dist = distance2(point, refPoints_[i]);

        // closer than previous closest point:
        if( dist < minDist )
        {
            minDist = dist;
            idxMinDist = i;
        }
    }

    // special case of last control point:
    if( idxMinDist == refPoints_.size() - 1 )
    {
        // simply map this to last interval:
        idxMinDist--;
    }

    // return index of interval of closest point:
    return idxMinDist;
}


/*!
 * Auxiliary function that maps a point in Cartesian coordinates onto an 
 * internal segment of the spline curve. This is accomplished by iteratively 
 * minimising the distance between the test point and a base point sampled from 
 * the spline curve using the (derivative free) method of Brent.
 *
 * \throws A logic error is thrown if the minimum can not be converged within 
 * a hardcoded number of 100 Brent iterations.
 */
gmx::RVec
SplineCurve3D::projectionInInterval(
        const gmx::RVec &point,
        const real &lo,
        const real &hi)
{
    // internal parameters:
    const boost::uintmax_t maxIter = 100;
    boost::uintmax_t iter = maxIter;

    // objective function binding:
    boost::function<real(real)> objFun;
    objFun = boost::bind(
            &SplineCurve3D::pointSqDist, 
            this, 
            point, 
            _1);
 
    // find minimum via Brent's method:
    int bits = std::numeric_limits<double>::digits;
    std::pair<double, double> result;
    result = boost::math::tools::brent_find_minima(
            objFun,
            lo, 
            hi,
            bits,
            iter);

    // make sure convergence has been reached:
    if( iter >= maxIter )
    {
        throw std::logic_error("Could not converge Brent iteration in "
                               "Cartesian to curvilinear mapping!");
    }

    // return curvilinear coordinates of point:
    // TODO: implement angular coordinate
    gmx::RVec curvPoint;
    curvPoint[SS] = result.first;
    curvPoint[RR] = result.second;
    return curvPoint;    
}


/*!
 * Auxiliary function that projects a point in Cartesian coordinates onto the
 * extrapolation range beyond its two endpoints. As the curve is known to be a
 * line in this range the projection is solved for analytically as a projection 
 * onto a ray. To achieve this, a ray is constructed by sampling two point 
 * from the spline curve (its endpoint and a second point \f$ ds \f$  beyond 
 * this endpoint. 
 */
gmx::RVec
SplineCurve3D::projectionInExtrapRange(
        const gmx::RVec &point,
        const real &ds)
{
    // return variable:
    gmx::RVec proj;

    // lower or upper extrapolation range?
    gmx::RVec extrapPointA;
    gmx::RVec extrapPointB;
    real arcLenOffset;
    real arcLenSign;
    if( ds < 0.0 )
    {
        // lower range:
        extrapPointA = ctrlPoints_.front();
        arcLenOffset = knots_.front();
        arcLenSign = -1.0;
        extrapPointB = this -> evaluate(arcLenOffset + ds, 0);
    }
    else if( ds > 0.0 )
    {
        // upper range:
        extrapPointA = ctrlPoints_.back();
        arcLenOffset = knots_.back();
        arcLenSign = 1.0;
        extrapPointB = this -> evaluate(arcLenOffset + ds, 0);
    }
    else
    {
        throw std::logic_error("Parameter dt may not be zero!");
    }

    // (non-normalised) direction vector for the extrapolating line:
    gmx::RVec lineDirVector;
    rvec_sub(extrapPointB, extrapPointA, lineDirVector);

    // vector between the test point and the ray's endpoint:
    gmx::RVec endpointVector;
    rvec_sub(point, extrapPointA, endpointVector);

    // is ray endpoint the closest point?
    // NOTE: in this case the angle between line direction vector and endpoint 
    // vector will be >= 90 degrees and the scalar product <= 0.0!
    real cosOfAngle = iprod(endpointVector, lineDirVector);
    if( cosOfAngle <= 0.0 )
    {
        proj[SS] = arcLenOffset;
        proj[RR] = distance2(point, extrapPointA);

        return proj;
    }

    // projection of test point position onto ray and base point:
    real b = cosOfAngle/iprod(lineDirVector, lineDirVector);
    gmx::RVec basePoint;
    svmul(b, lineDirVector, basePoint);
    rvec_add(basePoint, extrapPointA, basePoint);

    // curvilinear coordinates around this line:
    proj[SS] = arcLenOffset + arcLenSign*b;
    proj[RR] = distance2(point, basePoint);

    // return points in curvilinear coordinates:
    return proj;
}


/*!
 * Getter function for access to the spline curves control points.
 */
std::vector<gmx::RVec>
SplineCurve3D::ctrlPoints() const
{
    return ctrlPoints_;
}


/*!
 * Uses Newton-Cotes quadrature of curve speed to determine the length of the 
 * arc between two given parameter values. The specific quadrature rule applied
 * is Boole's rule.
 */
real
SplineCurve3D::arcLengthBoole(const real &lo, const real &hi)
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


/*!
 * Prepares a lookup table that associates an arc length with each knot, where
 * the first knot is assigned an arc length of zero.
 */
void
SplineCurve3D::prepareArcLengthTable()
{
    arcLengthTable_.resize(knots_.size());
    real segmentLength;
    for(unsigned int i = 0; i < knots_.size() - 1; i++)
    {
        // calculate length of current segment:
        segmentLength = arcLengthBoole(knots_[i], knots_[i+1]);

        // add to arc length table:
        arcLengthTable_[i + 1] = arcLengthTable_[i] + segmentLength;
    }

    // set flag:
    arcLengthTableAvailable_ = true;
}


/*!
 * Returns arc length value at the control points by simply removing the 
 * repeated knot values from the arc length lookup lookup table.
 */
std::vector<real>
SplineCurve3D::ctrlPointArcLength()
{
    // check availability:
    if( arcLengthTableAvailable_ == false )
    {
        prepareArcLengthTable();
    }

    // create copy of table and remove redundant elements:
    std::vector<real> arcLength(arcLengthTable_.begin() + degree_, 
                                arcLengthTable_.end() - degree_);

    return arcLength;
}


/*!
 * Returns arc length value at first control point.
 */
real
SplineCurve3D::frstPointArcLength()
{
    // check availability:
    if( arcLengthTableAvailable_ == false )
    {
        prepareArcLengthTable();
    }

    return arcLengthTable_.front();
}


/*!
 * Returns arc length value at last control point.
 */
real
SplineCurve3D::lastPointArcLength()
{
    // check availability:
    if( arcLengthTableAvailable_ == false )
    {
        prepareArcLengthTable();
    }

    return arcLengthTable_.back();
}


/*!
 * Returns the parameter value (in the current parameterisation, typically 
 * chord length) that corresponds to a given value of arc length. This is done
 * via bisection.
 */
real
SplineCurve3D::arcLengthToParam(real &arcLength)
{
    boost::uintmax_t maxIter = 100;
    real absTol = 0.01*std::sqrt(std::numeric_limits<real>::epsilon());

    // sanity check for arc length table:
    if( arcLengthTableAvailable_ != true )
    {
        prepareArcLengthTable();
    }

    // handle upper endpoint:
    if( arcLength == arcLengthTable_.back() )
    {
        return knots_.back();
    }

    // find appropriate interval:
    std::pair<std::vector<real>::iterator, std::vector<real>::iterator> bounds;
    bounds = std::equal_range(arcLengthTable_.begin(), arcLengthTable_.end(), arcLength);
    int idxLo = bounds.second - arcLengthTable_.begin() - 1;
    int idxHi = bounds.second - arcLengthTable_.begin();


    // handle query outside table range:
    // TODO: add case for query below lower bound and test!
    if( bounds.second == arcLengthTable_.end() )
    {
        return knots_.back() + arcLength - arcLengthTable_.back();
    }
    if( bounds.second == arcLengthTable_.begin() )
    {
        std::cerr<<"ERROR: arc length below table value range!"<<std::endl;
        std::abort();
    }

    // initialise bisection interval and lower limit:
    real tLo = knots_[idxLo];
    real tHi = knots_[idxHi];
    real tLi = tLo;
   
    // target arc length within this interval:
    real targetIntervalLength = arcLength - arcLengthTable_[idxLo];

    // termination condition binding:
    boost::function<bool(real, real)> termCond;
    termCond = boost::bind(&SplineCurve3D::arcLengthToParamTerm, 
                           this, 
                           _1, 
                           _2,
                           absTol);

    // objective function binding:
    boost::function<real(real)> objFun;
    objFun = boost::bind(&SplineCurve3D::arcLengthToParamObj, 
                         this, 
                         tLi, 
                         _1, 
                         targetIntervalLength);

        
    std::pair<real, real> result;
    result = boost::math::tools::toms748_solve(objFun,
                                               tLo, 
                                               tHi,
                                               termCond,
                                               maxIter);
    

    // best guess if middle of bracketing interval:
    return 0.5*(result.first + result.second);
}


/*!
 * Termination condition for re-parameterisation optimisation. 
 */
bool
SplineCurve3D::arcLengthToParamTerm(real lo, real hi, real tol)
{
    return std::abs(hi - lo) <= tol;
}


/*!
 * Objective function for re-parameterisation optimisation.
 */
real
SplineCurve3D::arcLengthToParamObj(real lo, real hi, real target)
{
    return arcLengthBoole(lo, hi) - target;
}


/*!
 * Computes the squared Euclidean distance between some point and the point 
 * on the spline curve with the given parameter value. For efficiency, the
 * square root is not drawn, hence the squared distance is returned.
 */
double
SplineCurve3D::pointSqDist(gmx::RVec point, double eval)
{
    // evaluate spline:
    gmx::RVec splPoint = evaluate(eval, 0);

    // return squared distance:
    return (splPoint[XX] - point[XX])*(splPoint[XX] - point[XX]) +
           (splPoint[YY] - point[YY])*(splPoint[YY] - point[YY]) +
           (splPoint[ZZ] - point[ZZ])*(splPoint[ZZ] - point[ZZ]);
}

