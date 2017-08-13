#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>

#include "path-finding/abstract_probe_path_finder.hpp"


/*
 *
 */
AbstractProbePathFinder::AbstractProbePathFinder(
        std::map<std::string, real> params,
        gmx::RVec initProbePos,
        t_pbc pbc,
        gmx::AnalysisNeighborhoodPositions porePos,
        std::vector<real> vdwRadii)
    : AbstractPathFinder(params)
    , vdwRadii_(vdwRadii)
    , initProbePos_(initProbePos)
    , crntProbePos_()
    , nbh_()
{

    // set parameters:
    if( params.find("pfProbeRadius") != params.end() )
    {
        probeRadius_ = params["pfProbeRadius"];
    }
    else
    {
        std::cerr<<"ERROR: No probe radius given!"<<std::endl;
        std::abort();
    }

    // find maximum vdw radius:
    maxVdwRadius_ = *std::max_element(vdwRadii.begin(), vdwRadii.end());
}


/*!
 * Sets parameters of the AnalysisNeighborhood object meinted by this class and
 * initialises an AnalysisneighborhoodSearch.
 */
void
AbstractProbePathFinder::prepareNeighborhoodSearch(
    t_pbc pbc,
    gmx::AnalysisNeighborhoodPositions porePos,
    real cutoff)
{
    // prepare analysis neighborhood:
    nbh_.setCutoff(cutoff);
    nbh_.setXYMode(false);
    nbh_.setMode(gmx::AnalysisNeighborhood::eSearchMode_Automatic);

    // initialise search:
    nbSearch_ = nbh_.initSearch(&pbc, porePos);
}


/*
 * Finds the minimal free distance, i.e. the shortest distance between the 
 * probe and the closest van-der-Waals surface.
 */
real
AbstractProbePathFinder::findMinimalFreeDistance(
        std::vector<real> optimSpacePos)
{
    // internal variables:
    real pairDist;              // distance between probe and pore atom
    real poreAtomVdwRadius;     // van-der-Waals radius of pore atom
    // TODO: using infinity here will cause a LAPACK error later in the code
    // IF the search cutoff is too small. Terminating the code in the case
    // may be a good idea overall, but better error handling is needed. Note 
    // that if an arbitraty value is chosen here the optimisation will still 
    // work evn with too small a cutoff, but that will cause problems later,
    // when the path points are interpolated: The very non-smooth spacing of
    // points will then lead to kinks in the spline!
    real minimalFreeDistance = std::numeric_limits<real>::infinity();            // radius of maximal non-overlapping sphere

    // convert point in optimisation space to point in configuration space:
    gmx::AnalysisNeighborhoodPositions probePos(optimToConfig(optimSpacePos).as_vec());

    // begin a pair search:
    gmx::AnalysisNeighborhoodPairSearch nbPairSearch = nbSearch_.startPairSearch(probePos);

    // loop over all pairs:
    gmx::AnalysisNeighborhoodPair pair;
    int i = 0;
    while( nbPairSearch.findNextPair(&pair) )
    {
        // get pair distance:
        // TODO: square root can be moved out of loop!
        pairDist = std::sqrt(pair.distance2());

        // get vdW radius of reference atom:
        poreAtomVdwRadius = vdwRadii_.at(pair.refIndex());

        // update void radius if necessary:
        // TODO: factor in probe radius!
        if( (pairDist - poreAtomVdwRadius - probeRadius_) < minimalFreeDistance )
        {
            minimalFreeDistance = pairDist - poreAtomVdwRadius;
        }
        i++;
    }
//    std::cout<<"i = "<<i<<std::endl;

    // return radius of maximal free sphere:
    return minimalFreeDistance; 
}

