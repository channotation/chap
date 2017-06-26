#ifndef MOLECULAR_PATH_HPP
#define MOLECULAR_PATH_HPP

#include <vector>
#include <map>


#include <gromacs/math/vec.h>
#include <gromacs/pbcutil/pbc.h>
#include <gromacs/utility/real.h>
#include <gromacs/selection/selection.h>

#include "geometry/spline_curve_1D.hpp"
#include "geometry/spline_curve_3D.hpp"


/*
 *
 */
class PathMappingParameters
{
    public:

        // parameters values:
        real nbhSearchCutoff_;
        real mapTol_;
        real extrapDist_;
        int numPathSamples_;
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

        // costructor and destructor:
        MolecularPath(std::vector<gmx::RVec> &pathPoints,
                      std::vector<real> &poreRadii);
        ~MolecularPath();

        // interface for mapping particles onto pathway:
        std::map<int, gmx::RVec> mapSelection(gmx::Selection mapSel,
                                              PathMappingParameters params,
                                              t_pbc *nbhSearchPbc);
        
        // check if points lie inside pore:
        std::map<int, bool> checkIfInside(
                const std::map<int, gmx::RVec> &mappedCoords,
                real margin);
        std::map<int, bool> checkIfInside(
                const std::map<int, gmx::RVec> &mappedCoords,
                real margin, 
                real sLo,
                real sHi);

        // access original points:
        std::vector<gmx::RVec> pathPoints();
        std::vector<real> pathRadii();

        // access properties of path:
        real length();
        std::pair<real, real> minRadius();
        real volume();
        real radius(real);
        real sLo();
        real sHi();

        // sample points from centreline:
        std::vector<real> sampleArcLength(int nPoints, real extrapDist);
        std::vector<gmx::RVec> samplePoints(int nPoints, real extrapDist);
        std::vector<gmx::RVec> samplePoints(std::vector<real> arcLengthSample);
        std::vector<gmx::RVec> sampleTangents(int nPoints, real extrapDist);
        std::vector<gmx::RVec> sampleTangents(std::vector<real> arcLengthSample);
        std::vector<gmx::RVec> sampleNormTangents(int nPoints, real extrapDist);
        std::vector<gmx::RVec> sampleNormTangents(std::vector<real> arcLengthSample);
        std::vector<gmx::RVec> sampleNormals(int nPoints, real extrapDist);
        std::vector<gmx::RVec> sampleNormals(std::vector<real> arcLengthSample);
        std::vector<real> sampleRadii(int nPoints, real extrapDist);
        std::vector<real> sampleRadii(std::vector<real> arcLengthSample);


    private:

        // mathematical constants:
        const real PI_ = std::acos(-1.0);

        // utilities for sampling functions:
        inline real sampleArcLenStep(int nPoints, real extrapDist); 

        // original path points and corresponding radii:
        std::vector<gmx::RVec> pathPoints_;
        std::vector<real> pathRadii_;

        // pore centre line and corresponding radius:
        SplineCurve3D centreLine_;
        SplineCurve1D poreRadius_;

        // properties of path:
        real openingLo_;
        real openingHi_;
        real length_;
};


#endif

