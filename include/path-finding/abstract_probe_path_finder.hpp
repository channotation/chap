#ifndef ABSTRACT_PROBE_PATH_FINDER
#define ABSTRACT_PROBE_PATH_FINDER

#include <vector>

#include <gromacs/trajectoryanalysis.h>
#include <gromacs/selection/nbsearch.h>

#include "path-finding/abstract_path_finder.hpp"
#include "path-finding/molecular_path.hpp"


/*
 * Abstract class that implements infrastructure used by all probe-based path
 * finding algorithms (such as the probe position).
 */
class AbstractProbePathFinder : public AbstractPathFinder
{
    public:

        // constructor:
        AbstractProbePathFinder(std::map<std::string, real> params,
                                gmx::RVec initProbePos,
                                std::vector<real> vdwRadii);


    protected:

        // auxiliary functions for setting neighborhood search parameters:
        void prepareNeighborhoodSearch(
                t_pbc *pbc,
                gmx::AnalysisNeighborhoodPositions porePos,
                real cutoff);


        int maxProbeSteps_;
        real probeStepLength_;
        real probeRadius_;
        real maxProbeRadius_;
        real nbhCutoff_;

        std::vector<real> vdwRadii_;
        real maxVdwRadius_;

        gmx::RVec initProbePos_;
        gmx::RVec crntProbePos_;

        t_pbc pbc_;
        gmx::AnalysisNeighborhood nbh_;
        gmx::AnalysisNeighborhoodSearch nbSearch_;
        
        real findMinimalFreeDistance(std::vector<real> optimSpacePos);

        // conversion between optimisation space and configruation space:
        virtual gmx::RVec optimToConfig(std::vector<real> optimSpacePos) = 0;
};


#endif

