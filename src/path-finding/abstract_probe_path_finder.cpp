#include <iostream>
#include <limits>

#include "path-finding/abstract_probe_path_finder.hpp"

/*
 * Constructor.
 *
 * TODO: nbSearch is currently passed as a raw pointer. This avoids copying, 
 * but is also rather not so elegent. Might be better to use a unique_ptr or 
 * shared_ptr? Maybe ask Gromacs mailing list about this?
 */
AbstractProbePathFinder::AbstractProbePathFinder(real probeStepLength,
                                                 real probeRadius,
                                                 real maxFreeDist,
                                                 int maxProbeSteps,
                                                 gmx::RVec &initProbePos,
                                                 std::vector<real> &vdwRadii,
                                                 gmx::AnalysisNeighborhoodSearch *nbSearch,
                                                 int saRandomSeed,
                                                 int saMaxCoolingIter,
                                                 int saNumCostSamples,
                                                 real saXi,
                                                 real saConvRelTol,
                                                 real saInitTemp,
                                                 real saCoolingFactor,
                                                 real saStepLengthFactor,
                                                 bool saUseAdaptiveCandGen)
    : AbstractPathFinder()
    , maxProbeSteps_(maxProbeSteps)
    , probeStepLength_(probeStepLength)
    , probeRadius_(probeRadius)
    , maxProbeRadius_(maxFreeDist)
    , vdwRadii_(vdwRadii)
    , initProbePos_(initProbePos)
    , crntProbePos_()
    , nbSearch_(nbSearch)
    , saRandomSeed_(saRandomSeed)
    , saMaxCoolingIter_(saMaxCoolingIter)
    , saNumCostSamples_(saNumCostSamples)
    , saXi_(saXi)
    , saConvRelTol_(saConvRelTol)
    , saInitTemp_(saInitTemp)
    , saCoolingFactor_(saCoolingFactor)
    , saStepLengthFactor_(saStepLengthFactor)
    , saUseAdaptiveCandGen_(saUseAdaptiveCandGen)
{

}

/*
 * Finds the minimal free distance, i.e. the shortest distance between the 
 * probe and the closest van-der-Waals surface.
 */
real
AbstractProbePathFinder::findMinimalFreeDistance(real *optimSpacePos)
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

    // prepare neighbourhood search:
    gmx::AnalysisNeighborhoodPairSearch nbPairSearch = nbSearch_ -> startPairSearch(probePos);
    gmx::AnalysisNeighborhoodPair pair;

    // loop over all pairs:
    while( nbPairSearch.findNextPair(&pair) )
    {
        // get pair distance:
        // TODO: square root can be moved out of loop!
        pairDist = std::sqrt(pair.distance2());

        // get vdW radius of reference atom:
        // TODO: factor in vdW radius!
        poreAtomVdwRadius = vdwRadii_.at(pair.refIndex());
        //poreAtomVdwRadius = 0.0;

        // update void radius if necessary:
        // TODO: factor in probe radius!
        if( (pairDist - poreAtomVdwRadius - probeRadius_) < minimalFreeDistance )
        {
            minimalFreeDistance = pairDist - poreAtomVdwRadius;
        }
    }

//    std::cout<<"minimalFreeDistance = "<<minimalFreeDistance<<std::endl;

    // return radius of maximal free sphere:
    return minimalFreeDistance; 
}
