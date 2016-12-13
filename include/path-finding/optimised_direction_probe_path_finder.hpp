#ifndef OPTIMISED_DIRECTION_PROBE_PATH_FINDER_HPP
#define OPTIMESED_DIRECTION_PROBE_PATH_FINDER_HPP

#include "path-finding/abstract_probe_path_finder.hpp"

/*
 *
 */
class OptimisedDirectionProbePathFinder : public AbstractProbePathFinder
{
    public:
        
        OptimisedDirectionProbePathFinder(real probeStepLength,
                                          real probeRadius,
                                          real maxFreeDist,
                                          int maxProbeSteps,
                                          gmx::RVec initProbePos,
                                          std::vector<real> vdwRadii,
                                          gmx::AnalysisNeighborhoodSearch *nbSearch,
                                          int saRandomSeed,
                                          int saMaxCoolingIter,
                                          int saNumCostSamples,
                                          real saXi,
                                          real saConvRelTol,
                                          real saInitTemp,
                                          real saCoolinGFactor,
                                          real saStepLengthFactor,
                                          bool saUseAdaptiveCandGen);

        void findPath();

    private:

        void optimiseInitialPos();
        void advanceAndOptimise(bool forward);

        gmx::RVec optimToConfig(real *optimSpacePos);
};

#endif

