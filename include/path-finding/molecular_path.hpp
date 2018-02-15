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


#ifndef MOLECULAR_PATH_HPP
#define MOLECULAR_PATH_HPP

#include <map>
#include <string>
#include <vector>

#include <gromacs/math/vec.h>
#include <gromacs/pbcutil/pbc.h>
#include <gromacs/utility/real.h>
#include <gromacs/selection/selection.h>

#include "external/rapidjson/document.h"

#include "geometry/spline_curve_1D.hpp"
#include "geometry/spline_curve_3D.hpp"


/*!
 * Enum for different methods for aligning molecular pathways between frames.
 */
enum ePathAlignmentMethod {ePathAlignmentMethodNone, 
                           ePathAlignmentMethodIpp};


/*! 
 * \brief Parameter container for MolecularPath path mapping functionality.
 *
 * Note that this is not currently used, but is retained in case parameters 
 * need to be reintroduced at a later time.
 */
class PathMappingParameters
{
    public:

};


/*!
 * \brief Representation of a molecular pathway.
 *
 * This class is used to describe a molecular pathway (e.g. an ion channel 
 * pore). It is typically created by a class derived from AbstractPathFinder
 * and permits access the pathways properties, such as its length(), volume(),
 * or minRadius(). A MolecularPathObjExporter can be used to generate a
 * mesh representing the pathways surface in the Wavefront Object format.
 *
 * Internally, the pathway is represented by a SplineCurve3D object, which 
 * describes a \f$ C^2 \f$ -continuous curve corresponding to the centre line
 * of the pathway. A further SplineCurve1D object is used to describe the 
 * pathway's radius along the centre line. Together, these splines provide a 
 * means of determining where in the pathway a particle is located using
 * mapSelection() and to decide whether a given particle lies inside the 
 * pathway or not using checkIfInside().
 *
 * The class also exposes several auxiliary functions such as samplePoints() or
 * sampleRadii() to provide access to the properties of the centre line curve
 * and radius spline directly.
 */
class MolecularPath
{
    public:

        // costructors and destructor:
        MolecularPath(
                std::vector<gmx::RVec> &pathPoints,
                std::vector<real> &poreRadii);
        MolecularPath(
                const rapidjson::Document &doc);
                /*
        MolecularPath(
                std::vector<real> poreRadiusKnots,
                std::vector<real> poreRadiusCtrlPoints,
                std::vector<real> centreLineKnots,
                std::vector<gmx::RVec> centreLineCtrlPoints);*/
        ~MolecularPath();

        // interface for mapping particles onto pathway:
        std::vector<gmx::RVec> mapPositions(
                const std::vector<gmx::RVec> &positions);
        std::map<int, gmx::RVec> mapSelection(
                const gmx::Selection &mapSel); 
        
        // check if points lie inside pore:
        std::map<int, bool> checkIfInside(
                const std::map<int, gmx::RVec> &mappedCoords,
                real margin);
        std::map<int, bool> checkIfInside(
                const std::map<int, gmx::RVec> &mappedCoords,
                real margin, 
                real sLo,
                real sHi);

        // centreline-mapped properties:
        void addScalarProperty(
                std::string name,
                SplineCurve1D property,
                bool divergent);
        std::map<std::string, std::pair<SplineCurve1D, bool>> scalarProperties() const;

        // access original points:
        std::vector<gmx::RVec> pathPoints();
        std::vector<real> pathRadii();

        // access internal spline curves:
        SplineCurve1D pathRadius();
        SplineCurve3D centreLine();

        // access aggregate properties of path:
        real length() const;
        std::pair<real, real> minRadius();
        real volume();
        real radius(real);
        real sLo();
        real sHi();

        // access properties of splines
        std::vector<real> poreRadiusKnots() const;
        std::vector<real> poreRadiusUniqueKnots() const;
        std::vector<real> poreRadiusCtrlPoints() const;
        std::vector<real> centreLineKnots() const;
        std::vector<real> centreLineUniqueKnots() const;
        std::vector<gmx::RVec> centreLineCtrlPoints() const;

        // sample points from centreline:
        std::vector<real> sampleArcLength(
                size_t nPoints, 
                real extrapDist) const;
        std::vector<gmx::RVec> samplePoints(
                size_t nPoints, 
                real extrapDist);
        std::vector<gmx::RVec> samplePoints(
                std::vector<real> arcLengthSample);
        std::vector<gmx::RVec> sampleTangents(
                size_t nPoints, real extrapDist);
        std::vector<gmx::RVec> sampleTangents(
                std::vector<real> arcLengthSample);
        std::vector<gmx::RVec> sampleNormTangents(
                size_t nPoints, 
                real extrapDist);
        std::vector<gmx::RVec> sampleNormTangents(
                std::vector<real> arcLengthSample);
        std::vector<gmx::RVec> sampleNormals(
                size_t nPoints, 
                real extrapDist);
        std::vector<gmx::RVec> sampleNormals(
                std::vector<real> arcLengthSample);
        std::vector<real> sampleRadii(
                size_t nPoints, 
                real extrapDist);
        std::vector<real> sampleRadii(
                std::vector<real> arcLengthSample);

        // change internal coordinate representation of path:
        void shift(const gmx::RVec &shift);
        

    private:

        // mathematical constants:
        const real PI_ = std::acos(-1.0);

        // utilities for sampling functions:
        inline real sampleArcLenStep(
                size_t nPoints, 
                real extrapDist) const; 

        // utilities for path mapping:
        inline gmx::RVec mapPosition(
                const gmx::RVec &cartCoord,
                const std::vector<real> &arcLenSample,
                const std::vector<gmx::RVec> &pathPointSample,
                const real mapTol);
        inline int numSamplePoints(const PathMappingParameters &params);

        // original path points and corresponding radii:
        std::vector<gmx::RVec> pathPoints_;
        std::vector<real> pathRadii_;

        // pore centre line and corresponding radius:
        SplineCurve3D centreLine_;
        SplineCurve1D poreRadius_;

        // properties mapped onto path:
        std::map<std::string, std::pair<SplineCurve1D, bool>> properties_;

        // properties of path:
        real openingLo_;
        real openingHi_;
        real length_;
};

#endif

