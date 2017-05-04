#include <iostream>

#include <gromacs/pbcutil/pbc.h>
#include <gromacs/selection/nbsearch.h>
#include <gromacs/selection/selection.h>

#include "geometry/cubic_spline_interp_3D.hpp"
#include "geometry/spline_curve_3D.hpp"

#include "path-finding/molecular_path.hpp"


/*
 * Constructor.
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

    // reparameterise centre line spline by arc length:
    centreLine_.arcLengthParam();

    // TODO: also interpolate radius appropriately!
}


/*
 * Destructor.
 */
MolecularPath::~MolecularPath()
{

}


/*
 * Maps all positions onto molecular pathway.
 */
std::map<int, gmx::RVec>
MolecularPath::mapSelection(gmx::Selection mapSel,
                            t_pbc *nbhSearchPbc)
{
    // create a set of reference positions on the pore centre line:
    // TODO: fine grain the reference point set
    gmx::AnalysisNeighborhoodPositions centreLinePos(pathPoints_);

    // prepare neighborhood search:
    real nbhSearchCutoff = 1.0;
    real mapTol = 1e-3;
    gmx::AnalysisNeighborhood nbh;
    nbh.setCutoff(nbhSearchCutoff);

    // create a neighborhood search with centre line points as reference:
    gmx::AnalysisNeighborhoodSearch nbhSearch = nbh.initSearch(nbhSearchPbc, centreLinePos);

    // build map of pathway mapped coordinates:
    std::map<int, gmx::RVec> mappedCoords;
    for(int i = 0; i < mapSel.posCount(); i++)
    {
        // find closest reference point on centre line:
        gmx::AnalysisNeighborhoodPair pair = nbhSearch.nearestPoint(mapSel.position(i));

        // check if reference point was found within cutoff distance:
        if( pair.isValid() == false )
        {
            continue;
        }

        // refine mapping:
        gmx::RVec mappedCoord = centreLine_.cartesianToCurvilinear(mapSel.position(i).x(),
                                                                   pair.refIndex(),
                                                                   mapTol);

        // TODO: check if this is within the local pore radius!

        // add to list of mapped coordinates:
        mappedCoords[mapSel.position(i).refId()] = mappedCoord;
    }

    return mappedCoords;
}


























