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

        AbstractProbePathFinder(std::map<std::string, real> params,
                                gmx::RVec initProbePos,
//                                gmx::AnalysisNeighborhoodSearch *nbSearch,
                                t_pbc pbc,
                                gmx::AnalysisNeighborhoodPositions porePos,
                                std::vector<real> vdwRadii);

    protected:



        int maxProbeSteps_;
        real probeStepLength_;
        real probeRadius_;
        real maxProbeRadius_;
        std::vector<real> vdwRadii_;

        gmx::RVec initProbePos_;
        gmx::RVec crntProbePos_;

        t_pbc pbc_;
        gmx::AnalysisNeighborhood nbh_;
        gmx::AnalysisNeighborhoodSearch nbSearch_;
        
        real findMinimalFreeDistance(std::vector<real> optimSpacePos);
        // TODO: need a second freeDistance function here that allows a different
        // implementation for the initial optimisation problem

        virtual gmx::RVec optimToConfig(std::vector<real> optimSpacePos) = 0;
};


#endif

