#ifndef INPLANE_OPTIMISED_PROBE_PATH_FINDER_HPP
#define INPLANE_OPTIMISED_PROBE_PATH_FINDER_HPP

#include <gromacs/trajectoryanalysis.h>

#include "path-finding/abstract_probe_path_finder.hpp"

/*
 *
 */
class InplaneOptimisedProbePathFinder : public AbstractProbePathFinder
{
    public:

        InplaneOptimisedProbePathFinder(real probeStepLength,
                                        gmx::RVec initProbePos,
                                        std::vector<real> vdwRadii,
                                        gmx::AnalysisNeighborhoodSearch nbSearch);

        void findPath();

    private:

        gmx::RVec chanDirVec_;
        gmx::RVec orthVecU_;
        gmx::RVec orthVecW_;

        void optimiseInitialPos();
        void advanceAndOptimise(bool forward);

        gmx::RVec optimToConfig(real *optimSpacePos);
};

#endif

