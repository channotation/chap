// CHAP - The Channel Annotation Package
// 
// Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
// Stephen J. Tucker
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


#ifndef INPLANE_OPTIMISED_PROBE_PATH_FINDER_HPP
#define INPLANE_OPTIMISED_PROBE_PATH_FINDER_HPP

#include <map>
#include <string>
#include <vector>

#include <gromacs/trajectoryanalysis.h>

#include "path-finding/abstract_probe_path_finder.hpp"


/*!
 * \brief Probe-based path-finder based on the HOLE algorithm.
 */
class InplaneOptimisedProbePathFinder : public AbstractProbePathFinder
{
    public:

        // constructor
        InplaneOptimisedProbePathFinder(std::map<std::string, real> params,
                                        gmx::RVec initProbePos,
                                        gmx::RVec chanDirVec,
                                        t_pbc *pbc,
                                        gmx::AnalysisNeighborhoodPositions porePos,
                                        std::vector<real> vdwRadii);

        // interface for setting parameters:
        void setParameters(const PathFindingParameters &params);

        // public interface for path finding:
        void findPath();

    private:

        gmx::AnalysisNeighborhoodPositions porePos_;
        t_pbc *pbc_;

        gmx::RVec chanDirVec_;
        gmx::RVec orthVecU_;
        gmx::RVec orthVecW_;

        void optimiseInitialPos();
        void advanceAndOptimise(bool forward);

        gmx::RVec optimToConfig(std::vector<real> optimSpacePos);
};

#endif

