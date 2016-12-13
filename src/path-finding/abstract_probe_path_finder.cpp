#include <iostream>

#include "path-finding/abstract_probe_path_finder.hpp"

/*
 *
 */
AbstractProbePathFinder::AbstractProbePathFinder(real probeStepLength,
                                                 gmx::RVec &initProbePos,
                                                 std::vector<real> &vdwRadii,
                                                 gmx::AnalysisNeighborhoodSearch nbSearch)
    : AbstractPathFinder()
    , probeStepLength_(probeStepLength)
    , probeRadius_(0.0)
    , vdwRadii_(vdwRadii)
    , initProbePos_(initProbePos)
    , crntProbePos_()
    , nbSearch_(nbSearch)
    , saRandomSeed_(15011991)
    , saMaxCoolingIter_(1e3)
    , saNumCostSamples_(50)
    , saXi_(3.0)
    , saConvRelTol_(1e-10)
    , saInitTemp_(10.0)
    , saCoolingFactor_(0.99)
    , saStepLengthFactor_(0.01)
    , saUseAdaptiveCandGen_(false)
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
    // TODO: infinity
    real minimalFreeDistance = 99;            // radius of maximal non-overlapping sphere

    // convert point in optimisation space to point in configuration space:
    gmx::AnalysisNeighborhoodPositions probePos(optimToConfig(optimSpacePos).as_vec());

    // prepare neighbourhood search:
    gmx::AnalysisNeighborhoodPairSearch nbPairSearch = nbSearch_.startPairSearch(probePos);
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
