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
                                gmx::RVec &initProbePos,
                                std::vector<real> &vdwRadii,
                                gmx::AnalysisNeighborhoodSearch nbSearch);

//    protected:

        real probeStepLength_;

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

        gmx::AnalysisNeighborhoodSearch nbSearch_;

        real findMinimalFreeDistance(real *optimSpacePos);

        virtual gmx::RVec optimToConfig(real *optimSpacePos) = 0;
};


#endif

