#ifndef ABSTRACT_PROBE_PATH_FINDER
#define ABSTRACT_PROBE_PATH_FINDER

#include <vector>

#include <gromacs/trajectoryanalysis.h>

#include "path-finding/abstract_path_finder.hpp"

/*
 * Abstract class that implements infrastructure used by all probe-based path
 * finding algorithms (such as the probe position).
 */
class AbstractProbePathFinder : public AbstractPathFinder
{
    public:

        AbstractProbePathFinder(real probeStepLength, 
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
                                bool saUseAdaptiveCandGen);

//    protected:

        int maxProbeSteps_;

        real probeStepLength_;
        real probeRadius_;
        real maxProbeRadius_;

        int saRandomSeed_;
        int saMaxCoolingIter_;
        int saNumCostSamples_;

        real saXi_;
        real saConvRelTol_;
        real saInitTemp_;
        real saCoolingFactor_;
        real saStepLengthFactor_;
        
        bool saUseAdaptiveCandGen_;

        std::vector<real> vdwRadii_;

        gmx::RVec initProbePos_;
        gmx::RVec crntProbePos_;

        gmx::AnalysisNeighborhoodSearch *nbSearch_;

        real findMinimalFreeDistance(real *optimSpacePos);
        // TODO: need a second freeDistance function here that allows a different
        // implementation for the initial optimisation problem

        virtual gmx::RVec optimToConfig(real *optimSpacePos) = 0;
};


#endif

