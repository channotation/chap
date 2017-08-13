#ifndef INPLANE_OPTIMISED_PROBE_PATH_FINDER_HPP
#define INPLANE_OPTIMISED_PROBE_PATH_FINDER_HPP

#include <map>
#include <string>
#include <vector>

#include <gromacs/trajectoryanalysis.h>

#include "path-finding/abstract_probe_path_finder.hpp"


/*
 *
 */
class InplaneOptimisedProbePathFinder : public AbstractProbePathFinder
{
    public:

        // constructor
        InplaneOptimisedProbePathFinder(std::map<std::string, real> params,
                                        gmx::RVec initProbePos,
                                        gmx::RVec chanDirVec,
                                        t_pbc pbc,
                                        gmx::AnalysisNeighborhoodPositions porePos,
                                        std::vector<real> vdwRadii);

        // interface for setting parameters:
        void setParameters(const PathFindingParameters &params);

        // public interface for path finding:
        void findPath();

    private:

        gmx::AnalysisNeighborhoodPositions porePos_;
        t_pbc pbc_;

        gmx::RVec chanDirVec_;
        gmx::RVec orthVecU_;
        gmx::RVec orthVecW_;

        void optimiseInitialPos();
        void advanceAndOptimise(bool forward);

        gmx::RVec optimToConfig(std::vector<real> optimSpacePos);
};

#endif

