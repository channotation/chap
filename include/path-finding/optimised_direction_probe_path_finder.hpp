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


#ifdef NOTDEFINED

#ifndef OPTIMISED_DIRECTION_PROBE_PATH_FINDER_HPP
#define OPTIMESED_DIRECTION_PROBE_PATH_FINDER_HPP

#include "path-finding/abstract_probe_path_finder.hpp"

/*!
 * \brief Not yet implemented.
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

        real probeStepLength_;

        gmx::RVec initProbePos_;
        gmx::RVec prevProbePos_;
        gmx::RVec crntProbePos_;

        matrix inverseRotationMatrix_;

        void optimiseInitialPos();
        void advanceAndOptimise(gmx::RVec initDirection);
        void updateInverseRotationMatrix(gmx::RVec direction);

        gmx::RVec optimToConfig(std::vector<real> optimSpacePos);
};

#endif

#endif
