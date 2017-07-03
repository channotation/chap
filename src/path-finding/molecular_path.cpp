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
 * between the given support points. Consqeuently, the input vectors must have
 * the same number of elements. This uses the arc length distance between 
 * centre line points as interpolation parameter.
 * 
 * Finally, the centre-line curve is reparameterised in terms of arc length
 * and the length of the path is computed as the arc length distance between
 * the first and last centre line points.
 */
MolecularPath::MolecularPath(std::vector<gmx::RVec> &pathPoints, 
                             std::vector<real> &pathRadii)
    : centreLine_()
    , poreRadius_()
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

    // reparameterise centre line spline by arc length:
    centreLine_.arcLengthParam();
}


/*!
 * Destructor.
 */
MolecularPath::~MolecularPath()
{

}


/*!
 * Function for mapping one point in Cartesian coordinates to spline 
 * coordinates used internally by both mapPosition() and mapSelection().
 *
 * The mapping procedes in two steps: First, the closest of a set of sampling
 * points is found by direct searching in order to find an initial point for
 * a minimisation procedure that iteratively minimises the distance between 
 * the Cartesian point and the spline curve using cartesianToCurvilinear().
 *
 * The set of sample points is provided externally as it can be the same for 
 * mutliple mapped points and resampling in each step would be unneccessarily
 * costly. An exception is thrown if the closest point in the initial step is 
 * one of the endpoints of the sample.
 *
 * \todo The first step is relatively slow if more than ca. 1000 sample points 
 * have to be checked. May be worth to implement a tree-based search here.
 */
inline
gmx::RVec
MolecularPath::mapPosition(const gmx::RVec &cartCoord,
                           const std::vector<real> &arcLenSample,
                           const std::vector<gmx::RVec> &pathPointSample,
                           const real mapTol)
{
    // find closest sample point:
    std::vector<real> distances;
    distances.reserve(pathPointSample.size());
    for(size_t j = 0; j < pathPointSample.size(); j++)
    {
        distances.push_back( distance2(cartCoord, pathPointSample[j]) );
    }
    int idxMinDist = std::min_element(distances.begin(), distances.end()) - distances.begin();

    // refine mapping by distance minimisation:
    gmx::RVec mappedCoord = centreLine_.cartesianToCurvilinear(
            cartCoord,
            arcLenSample[idxMinDist - 1],
            arcLenSample[idxMinDist + 1],
            mapTol);

    // check that all points have been mapped to the interior of the spline sample:
    if( idxMinDist == 0 || idxMinDist == (pathPointSample.size() - 1) )
    {
        throw std::runtime_error("Particle mapped onto endpoint samples.");
    }
   
    // return mapped coordinate:
    return mappedCoord;
}


/*!
 * Auxiliary function for calculating the number of sample points used in the
 * initial step of particle to path mapping. The number of sample points is 
 * determined according to
 *
 * \f[
 *      N = \left\lceil \frac{L + 2d}{\Delta s} \right\rceil
 * \f]
 *
 * where \f$ L \f$ is the length of the path, \f$ d \f$ is the extrapolation 
 * at either endpoint and \f$ \Delta s \f$ is the target step length. 
 */
inline
int
MolecularPath::numSamplePoints(const PathMappingParameters &params)
{
    return std::ceil( (2.0*params.extrapDist_ + length()) / params.sampleStep_ );
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
MolecularPath::mapPositions(const std::vector<gmx::RVec> &positions,
                            const PathMappingParameters &params)
{
    // prepare sample points on centre line with maximum arc length distance:
    int nPathSamples = numSamplePoints(params);
    std::vector<real> arcLenSample = sampleArcLength(
            nPathSamples, 
            params.extrapDist_);
    std::vector<gmx::RVec> pathPointSample = samplePoints(arcLenSample);

    // map all input positions onto centre line:
    std::vector<gmx::RVec> mappedPositions;
    mappedPositions.reserve(positions.size());
    for(auto it = positions.begin(); it != positions.end(); it++)
    {
        try
        {
            mappedPositions.push_back(mapPosition(*it,
                                                  arcLenSample,
                                                  pathPointSample,
                                                  params.mapTol_) );
        }
        catch( const std::runtime_error& e )
        {
            // TODO: adaptively increase extrapolation distance to handle 
            // failure where position is mapped onto sample endpoints, for now
            // fatal failure will have to do
            std::cerr<<"ERROR: "<<e.what()<<std::endl;
            std::abort();
        }        
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
MolecularPath::mapSelection(const gmx::Selection &mapSel,
                            const PathMappingParameters &params)
{
    // create a set of reference positions on the pore centre line:
    int nPathSamples = numSamplePoints(params);
    std::vector<real> arcLenSample = sampleArcLength(
            nPathSamples, 
            params.extrapDist_);
    const std::vector<gmx::RVec> pathPointSample = samplePoints(arcLenSample);

    // build map of pathway mapped coordinates:
    std::map<int, gmx::RVec> mappedCoords;
    for(size_t i = 0; i < mapSel.posCount(); i++)
    {
        try
        {
            mappedCoords[mapSel.position(i).refId()] = mapPosition(
                    mapSel.position(i).x(),
                    arcLenSample,
                    pathPointSample,
                    params.mapTol_);
        }
        catch( const std::runtime_error &e )
        {
            // TODO: adaptively increase extrapolation distance to handle cases
            // where some particles are mapped onto sample endpoints; for now 
            // a hard failure will have to suffice
            std::cerr<<"ERROR: "<<e.what()<<std::endl;
            std::abort();
        }
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
        isInside[it -> first] = (it -> second[1] < (poreRadius_(evalPoint, 0, eSplineEvalDeBoor)) + margin);
    }

    // return assessment:
    return isInside;
}


/*
 *
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
        // optain mapped s value:
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
 * Returns length of the pathway, defined as the the arc length distance 
 * between the first and last control point
 */
real
MolecularPath::length()
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
    return poreRadius_(s, 0, eSplineEvalDeBoor);
}


/*
 * TODO: update this with shifts of coordinate
 */
real
MolecularPath::sLo()
{
    return 0.0;
}


/*
 * TODO: update this with shifts of coordinate
 */
real
MolecularPath::sHi()
{
    return length_;
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
    real sMin;
    real sMax;
    if( itMin > r.begin() && itMin < r.end() )
    {
        sMin = s[idxMin - 1];
        sMax = s[idxMin + 1];
    }
    else if( itMin == r.begin() )
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
 *  of a path not being truely circular and (ii) the dependency of radius on
 *  centre line parameter not being \f$ \mathcal{O}(s^3) \f$ neccessarily. The 
 *  former effect is likely stronger so that the volume estimate should be 
 *  viewed as a lower bound.
 */
real
MolecularPath::volume()
{
    // initialise number of intervals larger than number of spline intervals:
    // (this ensures that function is polynomial on each interval)
    int numIntervals = (poreRadius_.nKnots()) - 2*poreRadius_.degree() - 1;

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
MolecularPath::sampleArcLength(int nPoints,
                               real extrapDist)
{
    // get spacing of points in arc length:
    real arcLenStep = sampleArcLenStep(nPoints, extrapDist);

    // evaluate spline to obtain sample points:
    std::vector<real> arcLengthSample;
    arcLengthSample.reserve(nPoints);
    for(int i = 0; i < nPoints; i++)
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
MolecularPath::samplePoints(int nPoints,
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
        points.push_back( centreLine_(arcLengthSample[i], 0, eSplineEvalDeBoor) );
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
MolecularPath::sampleTangents(int nPoints, real extrapDist)
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
MolecularPath::sampleNormTangents(int nPoints, real extrapDist)
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
MolecularPath::sampleNormals(int nPoints, real extrapDist)
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
MolecularPath::sampleNormals(std::vector<real> arcLengthSample)
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
MolecularPath::sampleRadii(int nPoints,
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
        radii.push_back( poreRadius_(evalPoint, 0, eSplineEvalDeBoor) );
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
        radii.push_back( poreRadius_(arcLengthSample[i], 0, eSplineEvalDeBoor) );
    }

    // return vector of points:
    return radii;
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
MolecularPath::sampleArcLenStep(int nPoints, real extrapDist)
{
    // get spacing of points in arc length:
    return (this -> length() + 2.0*extrapDist)/(nPoints - 1);
}

