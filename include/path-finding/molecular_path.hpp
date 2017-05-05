#ifndef MOLECULAR_PATH_HPP
#define MOLECULAR_PATH_HPP

#include <vector>
#include <map>

#include <gromacs/utility/real.h>                                               
#include <gromacs/math/vec.h> 

#include "geometry/spline_curve_1D.hpp"
#include "geometry/spline_curve_3D.hpp"


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
                                              t_pbc *nbhSearchPbc);

        // access properties of path:
        real length();

        // sample points on centreline:
        std::vector<real> sampleArcLength(int nPoints, real extrapDist);
        std::vector<gmx::RVec> samplePoints(int nPoints, real extrapDist);
        std::vector<gmx::RVec> samplePoints(std::vector<real> arcLengthSample);
        std::vector<real> sampleRadii(int nPoints, real extrapDist);
        std::vector<real> sampleRadii(std::vector<real> arcLengthSample);



    private:

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

