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


/*
 *
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

        // access original points:
        std::vector<gmx::RVec> pathPoints();
        std::vector<real> pathRadii();

        // access properties of path:
        real length();
        std::pair<real, real> minRadius();
        real volume();
        real radius(real);

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

