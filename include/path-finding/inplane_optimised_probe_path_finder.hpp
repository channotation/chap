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
//                                        gmx::AnalysisNeighborhoodSearch *nbSearch,
                                        t_pbc pbc,
                                        gmx::AnalysisNeighborhoodPositions porePos,
                                        std::vector<real> vdwRadii);

        void findPath();

    private:

        gmx::RVec chanDirVec_;
        gmx::RVec orthVecU_;
        gmx::RVec orthVecW_;

        void optimiseInitialPos();
        void advanceAndOptimise(bool forward);

        gmx::RVec optimToConfig(std::vector<real> optimSpacePos);
};

#endif

