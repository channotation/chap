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
#include <iostream>
#include <functional>
#include <limits>
#include <ctime>

#include <boost/math/tools/minima.hpp>

#include <gromacs/pbcutil/pbc.h>
#include <gromacs/selection/nbsearch.h>
#include <gromacs/selection/selection.h>

#include "geometry/cubic_spline_interp_1D.hpp"
#include "geometry/cubic_spline_interp_3D.hpp"
#include "geometry/spline_curve_1D.hpp"
#include "geometry/spline_curve_3D.hpp"

#include "path-finding/molecular_path.hpp"


/*!
 * Constructor to generate a MolecularPath object from a set of centre line 
 * points and corresponding radii.
 *
 * This uses CubicSplineInterp3D to perform cubic spline interpolation and 
 * build a \f$ C^2 \f$-continuous curve connecting all given centre line 
 * points. This uses the Euclidean distance between centre line points as 
 * interpolation parameter.
 *
 * Subsequently, CubicSplineInterp1D is used to smoothly interpolate the radius
 * between the given support points. Consequently, the input vectors must have
 * the same number of elements. This uses the arc length distance between 
 * centre line points as interpolation parameter.
 * 
 * Finally, the centre-line curve is re-parameterised in terms of arc length
 * and the length of the path is computed as the arc length distance between
 * the first and last centre line points.
 */
MolecularPath::MolecularPath(std::vector<gmx::RVec> &pathPoints, 
                             std::vector<real> &pathRadii)
{
    // assign internal containers for original path data:
    pathPoints_ = pathPoints;
    pathRadii_ = pathRadii;

    // construct centre line spline by interpolation of path points:
    CubicSplineInterp3D Interp3D;
    centreLine_ = Interp3D(pathPoints_, eSplineInterpBoundaryHermite);

    // get arc length at original control points:
    std::vector<real> arcLen = centreLine_.ctrlPointArcLength();
    openingLo_ = arcLen.front();
    openingHi_ = arcLen.back();
    length_ = std::abs(openingHi_ - openingLo_);

    // interpolate radius:
    CubicSplineInterp1D Interp1D;
    poreRadius_ =  Interp1D(arcLen, pathRadii_, eSplineInterpBoundaryHermite);

    // re-parameterise centre line spline by arc length:
    centreLine_.arcLengthParam();
}


/*!
 * Constructor for creating MolcularPath from a JSON document.
 */
MolecularPath::MolecularPath(
        const rapidjson::Document &doc)
    : centreLine_()
    , poreRadius_()
{
    // make sure document is valid object:
    if( !doc.IsObject() )
    {
        throw std::logic_error("JSON document passed to MolecularPath"
        "constructor is not a valid JSON object.");
    }

    // make sure radius spline data is given:
    if( !doc.HasMember("molPathRadiusSpline") )
    {
        throw std::logic_error("JSON document passed to MolecularPath"
        "constructor does not have required member molPathRadiusSpline.");
    }
    if( !doc["molPathRadiusSpline"].HasMember("knots") )
    {
        throw std::logic_error("Could not find knots in molPathRadiusSpline.");
    }
    if( !doc["molPathRadiusSpline"].HasMember("ctrl") )
    {
        throw std::logic_error("Could not find ctrl in molPathRadiusSpline.");
    }

    // make sure centre line spline data is given:
    if( !doc.HasMember("molPathCentreLineSpline") )
    {
        throw std::logic_error("JSON document passed to MolecularPath"
        "constructor does not have required member molPathCentreLineSpline.");
    }
    if( !doc["molPathCentreLineSpline"].HasMember("knots") )
    {
        throw std::logic_error("Could not find knots in molPathCentreLineSpline.");
    }
    if( !doc["molPathCentreLineSpline"].HasMember("ctrlX") )
    {
        throw std::logic_error("Could not find ctrlX in molPathCentreLineSpline.");
    }
    if( !doc["molPathCentreLineSpline"].HasMember("ctrlY") )
    {
        throw std::logic_error("Could not find ctrlY in molPathCentreLineSpline.");
    }
    if( !doc["molPathCentreLineSpline"].HasMember("ctrlZ") )
    {
        throw std::logic_error("Could not find ctrlZ in malPathCentreLineSpline.");
    }

    // make sure original points are present:
    if( !doc.HasMember("molPathOrigPoints") )
    {
        throw std::logic_error("JSON document passed to MolecularPath"
        "constructor does not have required member molPathOrigPoints.");
    }
    if( !doc["molPathOrigPoints"].HasMember("x") )
    {
        throw std::logic_error("Could not find x in molPathRadiusSpline.");
    }
    if( !doc["molPathOrigPoints"].HasMember("y") )
    {
        throw std::logic_error("Could not find y in molPathRadiusSpline.");
    }
    if( !doc["molPathOrigPoints"].HasMember("z") )
    {
        throw std::logic_error("Could not find z in molPathRadiusSpline.");
    }
    if( !doc["molPathOrigPoints"].HasMember("r") )
    {
        throw std::logic_error("Could not find r in molPathRadiusSpline.");
    }

    // extract original point and radius data from JSON:
    for(size_t i = 0; i < doc["molPathOrigPoints"]["r"].Size(); i++)
    {
         pathPoints_.push_back(
                gmx::RVec(doc["molPathOrigPoints"]["x"][i].GetDouble(), 
                          doc["molPathOrigPoints"]["y"][i].GetDouble(),
                          doc["molPathOrigPoints"]["z"][i].GetDouble()));
         pathRadii_.push_back(doc["molPathOrigPoints"]["r"][i].GetDouble());
    }

    // extract pore radius spline from data:
    std::vector<real> poreRadiusKnots;
    std::vector<real> poreRadiusCtrlPoints;
    for(size_t i = 0; i < doc["molPathRadiusSpline"]["knots"].Size(); i++)
    {
        poreRadiusKnots.push_back( 
                doc["molPathRadiusSpline"]["knots"][i].GetDouble() );
        poreRadiusCtrlPoints.push_back( 
                doc["molPathRadiusSpline"]["ctrl"][i].GetDouble() );
    } 

    // add duplicate knots at endpoints:
    int poreRadiusSplineDegree = 3; // TODO: should not be hardcoded
    poreRadiusKnots.insert(
            poreRadiusKnots.end(),
            poreRadiusSplineDegree - 1,
            poreRadiusKnots.back());
    poreRadiusKnots.insert(
            poreRadiusKnots.begin(),
            poreRadiusSplineDegree - 1,
            poreRadiusKnots.front());

    // create radius spline curve:
    poreRadius_ = SplineCurve1D(
            poreRadiusSplineDegree,
            poreRadiusKnots,
            poreRadiusCtrlPoints);

    // extract centre line spline from data:
    std::vector<real> centreLineKnots;
    std::vector<gmx::RVec> centreLineCtrlPoints;
    for(size_t i = 0; i < doc["molPathCentreLineSpline"]["knots"].Size(); i++)
    {
        centreLineKnots.push_back( 
                doc["molPathCentreLineSpline"]["knots"][i].GetDouble() );
        centreLineCtrlPoints.push_back( 
                gmx::RVec(doc["molPathCentreLineSpline"]["ctrlX"][i].GetDouble(),
                          doc["molPathCentreLineSpline"]["ctrlY"][i].GetDouble(),
                          doc["molPathCentreLineSpline"]["ctrlZ"][i].GetDouble()));
    } 

    // add duplicate knots at endpoints:
    int centreLineSplineDegree = 3; // TODO: should not be hardcoded
    centreLineKnots.insert(
            centreLineKnots.end(),
            centreLineSplineDegree - 1,
            centreLineKnots.back());
    centreLineKnots.insert(
            centreLineKnots.begin(),
            centreLineSplineDegree - 1,
            centreLineKnots.front());

    // create centre line spline curve:
    centreLine_ = SplineCurve3D(
            centreLineSplineDegree,
            centreLineKnots,
            centreLineCtrlPoints);

    // set position of openings and pore length:
    openingLo_ = poreRadiusKnots.front();
    openingHi_ = poreRadiusKnots.back();
    length_ = openingHi_ - openingLo_;

    // sanity check:
    if( openingLo_ > openingHi_ )
    {
        throw std::logic_error("Pore opening coordinates out of order.");
    }
}


/*!
 * Destructor.
 */
MolecularPath::~MolecularPath()
{

}


/*!
 * Function for mapping a set of Cartesian positions onto the centre line 
 * spline curve.
 *
 * The return value is a vector of points in spline coordinates ordered in the 
 * same way as the input vector. Internally, this uses mapPosition() for each 
 * input position.
 */
std::vector<gmx::RVec>
MolecularPath::mapPositions(const std::vector<gmx::RVec> &positions)
{
    // map all input positions onto centre line:
    std::vector<gmx::RVec> mappedPositions;
    mappedPositions.reserve(positions.size());
    for(auto pos : positions)
    {
        mappedPositions.push_back(centreLine_.cartesianToCurvilinear(pos));
    }
 
    // return mapped positions:
    return mappedPositions;
}


/*!
 * Maps all positions in a selection onto molecular pathway.
 *
 * This function does essentially the same as mapPositions, except that the 
 * input may be provided as a selection of particles and the output will be 
 * a map associating each refId in the selection with a set of curvilinear 
 * coordinates.
 */
std::map<int, gmx::RVec>
MolecularPath::mapSelection(const gmx::Selection &mapSel)
{
    // build map of pathway mapped coordinates:
    std::map<int, gmx::RVec> mappedCoords;
    for(int i = 0; i < mapSel.posCount(); i++)
    {
        unsigned int idx = mapSel.position(i).refId();
        mappedCoords[idx] = centreLine_.cartesianToCurvilinear(
                mapSel.position(i).x());
    }

    // return mapped coordinates:
    return mappedCoords;
}


/*!
 * Checks if points described by a set of mapped coordinates lie within the 
 * MolecularPath. 
 *
 * The input is taken to be a map of points in centre line coordinates, i.e.
 * a tuple of \f$ (s_i, \rho_i, \phi_i) \f$ values for the \f$ i \f$ -th point,
 * where \f$ s \f$ is the distance along the centre line, \f$ \rho \f$ is the 
 * (orthogonal) distance from the centre line, and \f$ \phi \f$ an angular 
 * coordinate which is ignored in this function. To obtain a set of mapped
 * points the mapSelection() method can be used.
 *
 * The output is a map of booleans indicating whether or not a point lies 
 * inside the MolecularPath. The integer key is the same as in the input map
 * and will usually correspond to the ID of a particle.
 *
 * To test whether a point lies within the MolecularPath, its distance from 
 * the centre line is compared to to path radius at this point plus a margin,
 * i.e.
 *
 * \f[
 *      \gamma_i = \left\{ 
 *                 \begin{array}{ll}
 *                     1 & \text{ if } \rho_i < R(s_i) + m \\
 *                     0 & \text{ otherwise }
 *                 \end{array}
 *                 \right.
 * \f]
 *
 * with \f$ \gamma_i \f$ being a binary indicator function.
 */
std::map<int, bool>
MolecularPath::checkIfInside(const std::map<int, gmx::RVec> &mappedCoords,
                             real margin)
{
    // create map for check results:
    std::map<int, bool> isInside;

    for(auto it = mappedCoords.begin(); it != mappedCoords.end(); it++)
    {
        real evalPoint = it -> second[0];
        real thres = (poreRadius_.evaluate(evalPoint, 0)) + margin;
        // threshold needs to be squared here because radial coordinate is!
        isInside[it -> first] = (it -> second[1] < thres*thres);
    }

    // return assessment:
    return isInside;
}


/*!
 * Checks if a given set of coordinates lies within the pathway.
 */
std::map<int, bool>
MolecularPath::checkIfInside(const std::map<int, gmx::RVec> &mappedCoords,
                             real margin,
                             real sLo,
                             real sHi)
{
    // first make decision based on margin:
    std::map<int, bool> isInside = checkIfInside(mappedCoords, margin);

    // now erase all elements that do not fall in given range:
    for(auto it = isInside.begin(); it != isInside.end(); it++)
    {
        // obtain mapped s value:
        real s = mappedCoords.at(it->first)[0];

        // does it fall in target interval?
        if( s < sLo || s > sHi )
        {
            it -> second = false;
        }
    }

    // return resulting map:
    return isInside;
}


/*!
 * Adds a scalar property to the MolecularPath. Note that property names must
 * be unique and already existing properties will be overwritten.
 */
void
MolecularPath::addScalarProperty(
        std::string name,
        SplineCurve1D property,
        bool divergent)
{
    properties_[name] = std::pair<SplineCurve1D, bool>(property, divergent);
}


/*!
 * Returns map of scalar properties associated with the MolecularPath.
 */
std::map<std::string, std::pair<SplineCurve1D, bool>>
MolecularPath::scalarProperties() const
{
    return properties_;
}


/*! 
 * Simple getter function for access to original path points used to construct
 * the path.
 */
std::vector<gmx::RVec>
MolecularPath::pathPoints()
{
    return pathPoints_;
}


/*!
 * Simple getter function for access to original radii used to construct the 
 * path.
 */
std::vector<real>
MolecularPath::pathRadii()
{
    return pathRadii_;
}


/*!
 * Returns a copy of the internal pore radius spline.
 */
SplineCurve1D
MolecularPath::pathRadius()
{
    return poreRadius_;
}


/*!
 * Returns a copy of the internal centre line spline.
 */
SplineCurve3D
MolecularPath::centreLine()
{
    return centreLine_;
}


/*! 
 * Returns length of the pathway, defined as the the arc length distance 
 * between the first and last control point
 */
real
MolecularPath::length() const
{
    return length_;
}


/*!
 * Returns the pore radius \f$ R(s) \f$ at a given value of the centre line's
 * spline parameter \f$ s \f$.
 */
real
MolecularPath::radius(real s)
{
    return poreRadius_.evaluate(s, 0);
}


/*!
 * Returns coordinates of the lower opening of the pore.
 */
real
MolecularPath::sLo()
{
    return openingLo_;
}


/*!
 * Returns the coordinates of the upper opening of the pore.
 */
real
MolecularPath::sHi()
{
    return openingHi_;
}


/*!
 * Getter method for access to the radius spline's knot vector. This returns 
 * the complete knot vector including duplicate points at the ends.
 */
std::vector<real>
MolecularPath::poreRadiusKnots() const
{
    return poreRadius_.knotVector();
}


/*!
 * Getter method for access to the radius spline's knot vector. This does strip
 * the knot vector of repeated knots at end points so that the resulting vector 
 * has as many elements as the vector of control points.
 */
std::vector<real>
MolecularPath::poreRadiusUniqueKnots() const
{           
    std::vector<real> allKnots = poreRadius_.knotVector();
    std::vector<real> uniqueKnots(
        allKnots.begin() + poreRadius_.degree() - 1, 
        allKnots.end() - poreRadius_.degree() + 1);
    return uniqueKnots;
}


/*!
 * Getter method for access to the radius spline's control points.
 */
std::vector<real>
MolecularPath::poreRadiusCtrlPoints() const
{
    return poreRadius_.ctrlPoints();
}


/*!
 * Getter method for access to the centre line spline's knot vector. 
 * This returns  the complete knot vector including duplicate points at the 
 * ends.
 */
std::vector<real>
MolecularPath::centreLineKnots() const
{
    return centreLine_.knotVector();
}


/*!
 * Getter method for access to the centre line spline's knot vector. This does 
 * strip the knot vector of repeated knots at end points so that the resulting 
 * vector has as many elements as the vector of control points.
 */
std::vector<real>
MolecularPath::centreLineUniqueKnots() const
{
    std::vector<real> allKnots = centreLine_.knotVector();
    std::vector<real> uniqueKnots(
            allKnots.begin() + centreLine_.degree() - 1,
            allKnots.end() - centreLine_.degree() + 1);
    return uniqueKnots;
}


/*!
 * Getter method for access to the centre line spline's control points.
 */
std::vector<gmx::RVec>
MolecularPath::centreLineCtrlPoints() const
{
    return centreLine_.ctrlPoints();
}



/*!
 * Finds the minimum radius of the path and the location along the centre line
 * (in the current parameterisation) of this minimum. 
 *
 * This is achieved by first sampling a set of trial points no more than 0.1
 * nm apart along the centreline and evaluating the path radius at all of these
 * points. The position of the minimum radius sample point is then used to 
 * create input to Brent's minimisation algorithm. The initial bracketing 
 * interval for Brent's algorithm is taken as one sampling point above and 
 * below the minimum radius sampling point, if this point lies somewhere in 
 * between the the centre line's endpoints. If the minimum radius sampling
 * point falls on either endpoint of the centeline, this point itself is taken
 * as one of the bracketing interval's limits.
 *
 * Note that Brent's algorithm is currently limited to a hardcoded limit of
 * 100 iterations.
 */
std::pair<real, real>
MolecularPath::minRadius()
{
    // internal parameters:
    real maxSampleDist = 0.1;
    boost::uintmax_t maxIter = 100;

    // draw radius samples along path:
    int nSamples = std::ceil(length_/maxSampleDist);
    std::vector<real> s = sampleArcLength(nSamples, 0.0);
    std::vector<real> r = sampleRadii(s);
   
    // find smallest sample radius:
    auto itMin = std::min_element(r.begin(), r.end());
    int idxMin = std::distance(r.begin(), itMin);

    // determine bracketing interval:
    real sMin = s[idxMin - 1];
    real sMax = s[idxMin + 1];
    if( itMin == r.begin() )
    {
        sMin = s[idxMin];
        sMax = s[idxMin + 1];
    }
    else if( itMin == r.end() )
    {
        sMin = s[idxMin - 1];
        sMax = s[idxMin];
    }

    // return minimum and arg min:    
    return boost::math::tools::brent_find_minima(
            std::bind(&MolecularPath::radius, this, std::placeholders::_1), 
            sMin, 
            sMax, 
            std::numeric_limits<real>::digits,
            maxIter);
}


/*!
 *  Returns the volume of the path defined as the volume of the spline tube
 *  between the upper and lower opening of the path, i.e.
 *
 *  \f[
 *
 *      V = \pi \int_{s_0}^{s_1} \left( R(s) \right)^2 ds
 *
 *  \f]
 *
 *  where \f$ R(s) \f$ denotes the radius at a given point along the spline.
 *
 *  This integral is solved numerically by applying a sixth order Newton-Cotes
 *  scheme (Weddle's rule) in each interval between two subsequent knots. Since
 *  the path radius is \f$ \mathcal{O}(s^3) \f$ by construction, the integrand 
 *  is guaranteed to be of order \f$ \mathcal{O}(s^6) \f$ so the integration
 *  is exact.
 *
 *  The volume is nonetheless an estimate due to (i) the cross-sectional area
 *  of a path not being truly circular and (ii) the dependency of radius on
 *  centre line parameter not being \f$ \mathcal{O}(s^3) \f$ necessarily. The 
 *  former effect is likely stronger so that the volume estimate should be 
 *  viewed as a lower bound.
 */
real
MolecularPath::volume()
{
    // initialise number of intervals larger than number of spline intervals:
    // (this ensures that function is polynomial on each interval)
    size_t numIntervals = (poreRadius_.nKnots()) - 2*poreRadius_.degree() - 1;

    // integration interval:
    real h = length_/numIntervals;

    // sample path radii:
    int numPoints = 6*numIntervals + 1;
    std::vector<real> radii = sampleRadii(numPoints, 0.0);
    std::vector<real> s = sampleArcLength(numPoints, 0.0);
  
    // calculate squared radii:
    std::vector<real> sqRadii;
    sqRadii.reserve(radii.size());
    for(size_t i = 0; i < radii.size(); i++)
    {
        sqRadii.push_back(radii[i]*radii[i]);
    }

    // evaluate integral using fine grained support points:
    real integral = 0.0;
    for(size_t i = 0; i < numIntervals; i++)
    {
        // index to vector of radii:
        int idx = 6*i;

        // integrate this interval using Weddle's Rule:
        integral +=  41.0*sqRadii[idx] + 
                    216.0*sqRadii[idx+1] +
                     27.0*sqRadii[idx+2] +
                    272.0*sqRadii[idx+3] + 
                     27.0*sqRadii[idx+4] +
                    216.0*sqRadii[idx+5] +
                     41.0*sqRadii[idx+6];
    }
    integral *= h/840.0;

    // multiply in constant factor:
    integral *= PI_;

    // return value of integral:
    return integral;
}


/*!
 * Returns a vector of equally spaced arc length points that extends a 
 * specified distance beyond the openings of the pore.
 */
std::vector<real>
MolecularPath::sampleArcLength(size_t nPoints,
                               real extrapDist) const
{
    // get spacing of points in arc length:
    real arcLenStep = sampleArcLenStep(nPoints, extrapDist);

    // evaluate spline to obtain sample points:
    std::vector<real> arcLengthSample;
    arcLengthSample.reserve(nPoints);
    for(size_t i = 0; i < nPoints; i++)
    {
        // calculate evaluation point:
        arcLengthSample.push_back( openingLo_ - extrapDist + i*arcLenStep );  
    }

    // return vector of points:
    return arcLengthSample;
}



/*!
 * Returns a vector of point on the molecular path's centre line. The points 
 * will be equally spaced in arc length and sampling will extend beyond the 
 * pore openings for the specified distance. All points are given in
 * Cartesian coordinates.
 */
std::vector<gmx::RVec>
MolecularPath::samplePoints(size_t nPoints,
                            real extrapDist)
{
    // sample equidistant arc length values:
    std::vector<real> arcLengthSteps = sampleArcLength(nPoints, extrapDist);

    // return vector of points at these values:
    return samplePoints(arcLengthSteps);
}


/*! 
 * Returns a vector of points on the path's centre line at the given arc length
 * parameter values. All points are given in Cartesian coordinates.
 */
std::vector<gmx::RVec>
MolecularPath::samplePoints(std::vector<real> arcLengthSample)
{
    // evaluate spline to obtain sample points:
    std::vector<gmx::RVec> points;
    points.reserve(arcLengthSample.size());
    for(size_t i = 0; i < arcLengthSample.size(); i++)
    {
        // evaluate spline at this point:
        points.push_back( centreLine_.evaluate(arcLengthSample[i], 0) );
    }

    // return vector of points:
    return points;
}


/*!
 * Returns vector of \p nPoints tangents to the centre line. The samples are 
 * taken from equidistant points along the spline, extending \p extrapDist into 
 * the extrapolation range on either side. 
 */
std::vector<gmx::RVec>
MolecularPath::sampleTangents(size_t nPoints, real extrapDist)
{
    // sample equidistant arc length values:
    std::vector<real> arcLengthSteps = sampleArcLength(nPoints, extrapDist);

    // return vector of tangents at these values:
    return sampleTangents(arcLengthSteps);
}


/*!
 * Returns vector of tangents to the centre line. These are calculated at the 
 * evaluation points given on \p arcLengthSample.
 */
std::vector<gmx::RVec>
MolecularPath::sampleTangents(std::vector<real> arcLengthSample)
{    
    // evaluate spline to obtain sample points:
    std::vector<gmx::RVec> tangents;
    tangents.reserve(arcLengthSample.size());
    for(size_t i = 0; i < arcLengthSample.size(); i++)
    {
        // evaluate spline at this point:
        tangents.push_back( centreLine_.tangentVec(arcLengthSample[i]) );
    }

    // return vector of points:
    return tangents;
}


/*!
 * Returns vector of \p nPoints tangents to the centre line. The samples are 
 * taken from equidistant points along the spline, extending \p extrapDist into 
 * the extrapolation range on either side. All tangents are explicitly 
 * normalised in this function.
 */

std::vector<gmx::RVec>
MolecularPath::sampleNormTangents(size_t nPoints, real extrapDist)
{
    // sample equidistant arc length values:
    std::vector<real> arcLengthSteps = sampleArcLength(nPoints, extrapDist);

    // return vector of tangents at these values:
    std::vector<gmx::RVec> tangents = sampleTangents(arcLengthSteps);

    // normalise all tangent vectors:
    std::vector<gmx::RVec>::iterator it;
    for(it = tangents.begin(); it != tangents.end(); it++)
    {
        unitv(*it, *it);
    }

    return tangents;
}

/*!
 * Returns vector of tangents to the centre line. These are calculate at the 
 * evaluation points given on \p arcLengthSample. All tangents are explicitly
 * normalised in this function.
 */
std::vector<gmx::RVec>
MolecularPath::sampleNormTangents(std::vector<real> arcLengthSample)
{    
    // evaluate spline to obtain sample points:
    std::vector<gmx::RVec> tangents;
    tangents.reserve(arcLengthSample.size());
    for(size_t i = 0; i < arcLengthSample.size(); i++)
    {
        // evaluate spline at this point:
        tangents.push_back( centreLine_.tangentVec(arcLengthSample[i]) );

        // normalise tangent vector:
        unitv(tangents.back(), tangents.back());
    }

    // return vector of points:
    return tangents;
}


/*!
 * \todo This needs to be implemented.
 */
std::vector<gmx::RVec>
MolecularPath::sampleNormals(size_t nPoints, real extrapDist)
{
    // sample equidistant arc length values:
    std::vector<real> arcLengthSteps = sampleArcLength(nPoints, extrapDist);

    // return vector of normals at these values:
    return sampleNormals(arcLengthSteps);
}


/*!
 * \todo This needs to be implemented.
 */
std::vector<gmx::RVec>
MolecularPath::sampleNormals(std::vector<real> /* arcLengthSample */)
{
    // evaluate spline to obtain sample points:
    std::vector<gmx::RVec> normals;

    // return vector of points:
    return normals;
}


/*!
 * Returns a vector of radius values at equally spaced points a long the path.
 * Sampling extends the specified distance beyond the openings of the pore.
 */
std::vector<real>
MolecularPath::sampleRadii(size_t nPoints,
                           real extrapDist)
{
    // get spacing of points in arc length:
    real arcLenStep = sampleArcLenStep(nPoints, extrapDist);

    // evaluate spline to obtain sample points:
    std::vector<real> radii;
    for(size_t i = 0; i < nPoints; i++)
    {
        // calculate evaluation point:
        real evalPoint = openingLo_ - extrapDist + i*arcLenStep;  

        // evaluate spline at this point:
        radii.push_back( poreRadius_.evaluate(evalPoint, 0) );
    }

    // return vector of points:
    return radii;
}


/*!
 * Returns a vector of radius values at the given arc length parameter values.
 */
std::vector<real>
MolecularPath::sampleRadii(std::vector<real> arcLengthSample)
{
    // evaluate spline to obtain sample points:
    std::vector<real> radii;
    for(size_t i = 0; i < arcLengthSample.size(); i++)
    {
        // evaluate spline at this point:
        radii.push_back( poreRadius_.evaluate(arcLengthSample[i], 0) );
    }

    // return vector of points:
    return radii;
}


/*!
 * Shift the s-coordinate by the given number.
 */
void
MolecularPath::shift(const gmx::RVec &shift)
{
    // shift both the centre line and the radius spline:
    centreLine_.shift(shift);  
    poreRadius_.shift(shift);

    // adjust convenience variables defined in MolecularPath itself:
    openingLo_ -= shift[SS];
    openingHi_ -= shift[SS];
}


/*!
 * Auxiliary function that computes the step length along the arc for sampling
 * a given property at \f$ N \f$ points and reaching \f$ d \f$ into the 
 * extrapolation region at either side of the pore. The step length is computed
 * as
 *
 * \f[
 *      \Delta s = \frac{L + 2 d}{N - 1}
 * \f]
 *
 * where \f$ L \f$ denotes the length of the path between its first and last
 * control points.
 */
real
MolecularPath::sampleArcLenStep(size_t nPoints, real extrapDist) const
{
    // get spacing of points in arc length:
    return (this -> length() + 2.0*extrapDist)/(nPoints - 1);
}

