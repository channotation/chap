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


#include <iostream>
#include <limits>

#include <gromacs/math/vec.h>

#include "optim/simulated_annealing_module.hpp"
#include "optim/nelder_mead_module.hpp"

#include "path-finding/inplane_optimised_probe_path_finder.hpp"


/*!
 * Constructor.
 */
InplaneOptimisedProbePathFinder::InplaneOptimisedProbePathFinder(
        std::map<std::string, real> params,
        gmx::RVec initProbePos,
        gmx::RVec chanDirVec,
        t_pbc *pbc,
        gmx::AnalysisNeighborhoodPositions porePos,
        std::vector<real> vdwRadii)
    : AbstractProbePathFinder(params, initProbePos, vdwRadii)
    , porePos_(porePos)
    , pbc_(pbc)
    , chanDirVec_(chanDirVec)
    , orthVecU_(0.0, 0.0, 0.0)
    , orthVecW_(0.0, 0.0, 0.0)
{
    // tolerance threshold for norm of vector (which should be unit vectors):
    real nonZeroTol = std::numeric_limits<real>::epsilon();
    if( norm(chanDirVec_) < nonZeroTol )
    {
        throw std::runtime_error("Channel direction vector has norm close to "
                                 "zero. Please provide a finite-length channel "
                                 "direction vector with -pf-chan-dir-vec.");
    }

    // normalise channel direction vector:
    unitv(chanDirVec_, chanDirVec_);

    // generate first orthogonal vector:
    orthVecU_ = gmx::RVec(-chanDirVec_[YY], chanDirVec_[XX], 0.0);

    // make sure it is non null vector:
    if( norm(orthVecU_) < nonZeroTol )
    {
        // try different permutation:
        orthVecU_ = gmx::RVec(-chanDirVec_[ZZ], 0.0, chanDirVec_[XX]);

        // make sure it is not null vector:
        if( norm(orthVecU_) < nonZeroTol )
        {
            // try different permutation:
            orthVecU_ = gmx::RVec(0.0, -chanDirVec_[2], chanDirVec_[1]);

            // make sure it is not null vector:
            if( norm(orthVecU_) < nonZeroTol )
            {
                throw std::logic_error("Inplane optimised probe path finder "
                                       "could not generate an orthogonal "
                                       "vector.");
            }
        }
    }

    // normalise first orthogonal vector:
    unitv(orthVecU_, orthVecU_);

    // generate second orthogonal vector:
    // (this will already be normalised)
    cprod(chanDirVec_, orthVecU_, orthVecW_);

    // make sure vectors are mutually orthogonal:
    // (note that the vectors are not strictly required to be orthogonal, they 
    // may just not be colinear, so small numeric deviations should not be a
    // problem here)
    real orthTol = std::numeric_limits<real>::epsilon();
    if( iprod(chanDirVec_, orthVecU_) > orthTol )
    {
        throw std::logic_error("Vectors v and u in inplane optimised probe "
                               "path finder are not orthogonal.");
    }
    if( iprod(chanDirVec_, orthVecW_) > orthTol )
    {
        throw std::logic_error("Vectors v and w in inplane optimised probe "
                               "path finder are not orthogonal.");
    }
    if( iprod(orthVecU_, orthVecW_) > orthTol )
    {
        throw std::logic_error("Vectors u and w in inplane optimised probe "
                               "path finder are not orthogonal.");
    }
}


/*!
 * Set parameters for path-finding.
 */
void
InplaneOptimisedProbePathFinder::setParameters(
        const PathFindingParameters &params)
{
    // set parameters:
    probeStepLength_ = params.probeStepLength();
    maxProbeRadius_ = params.maxProbeRadius();
    maxProbeSteps_ = params.maxProbeSteps();

    // has cutoff been set by user:
    if( params.nbhCutoffIsSet() )
    {
        // user given cutoff:
        nbhCutoff_ = params.nbhCutoff();
    }
    else
    {
        // calculate cutoff automatically:
        real safetyMargin = std::sqrt(std::numeric_limits<real>::epsilon());
        nbhCutoff_ = params.maxProbeRadius() + maxVdwRadius_ + safetyMargin;
    }

    // set flag to true:
    parametersSet_ = true;
}


/*!
 * Execute path-finding algorithm.
 */
void
InplaneOptimisedProbePathFinder::findPath()
{
    // sanity check:
    if( !parametersSet_ )
    {
        throw std::logic_error("Path finding parameters have not been set.");
    }

    // prepare neighborhood search:
    prepareNeighborhoodSearch(
            pbc_,
            porePos_,
            nbhCutoff_);

    // optimise initial position:
    optimiseInitialPos();
    
    // advance forward:
    advanceAndOptimise(true);

    // revert array:
    std::reverse(path_.begin(), path_.end());
    std::reverse(radii_.begin(), radii_.end());
    
    // advance backward:
    advanceAndOptimise(false);
}


/*!
 * Optimise initial position of probe.
 */
void
InplaneOptimisedProbePathFinder::optimiseInitialPos()
{
    // set current probe position to initial probe position: 
    crntProbePos_ = initProbePos_;

    // initial state in optimisation space is always null vector:
    std::vector<real> initState = {0.0, 0.0};
    
    // cost function is minimal free distance function:
    ObjectiveFunction objFun;
    objFun = std::bind(&InplaneOptimisedProbePathFinder::findMinimalFreeDistance, 
                       this, std::placeholders::_1);

    // optimise in plane through simulated annealing:
    SimulatedAnnealingModule sam;
    sam.setObjFun(objFun);
    sam.setParams(params_);
    sam.setInitGuess(initState);
    sam.optimise();

    // refine with Nelder-Mead optimisation:
    NelderMeadModule nmm;
    nmm.setObjFun(objFun);
    nmm.setParams(params_);
    nmm.setInitGuess(sam.getOptimPoint().first);
    nmm.optimise();
       
    // set initial position to its optimal value:
    initProbePos_ = optimToConfig(nmm.getOptimPoint().first);

    // handle situation where cutoff radius was too small:
    // (or otherwise no particle was found within cutoff radius)
    if( std::isinf( nmm.getOptimPoint().second ) )
    {
        throw std::runtime_error("Pore radius at initial probe position is "
                                 "infinite. Consider increasing the maximum "
                                 "pore radius with -pf-max-free-dist or set "
                                 "an appropriate cutoff for neighbourhood "
                                 "searches explicitly with -pf-cutoff.");
    }

    // add path support point and associated radius to container:
    path_.push_back(initProbePos_);
    radii_.push_back(nmm.getOptimPoint().second);   
}


/*!
 * Optimise probe position in subsequent parallel planes.
 */
void
InplaneOptimisedProbePathFinder::advanceAndOptimise(bool forward)
{
    // set previous position to initial point:
    crntProbePos_ = initProbePos_;

    // set up direction vector for forward/backward marching:
    gmx::RVec direction(chanDirVec_);
    if( !forward )
    {
        direction[XX] = -direction[XX];
        direction[YY] = -direction[YY];
        direction[ZZ] = -direction[ZZ];
    }

    // initial state in optimisation space is always null vector:
    std::vector<real> initState = {0.0, 0.0};

    // cost function is minimal free distance function:
    ObjectiveFunction objFun;
    objFun = std::bind(&InplaneOptimisedProbePathFinder::findMinimalFreeDistance, 
                       this, std::placeholders::_1);


    // advance probe in direction of (inverse) channel direction vector:
    int numProbeSteps = 0;
    while(true)
    {
        // advance probe position to next plane:
        crntProbePos_[XX] = crntProbePos_[XX] + probeStepLength_*direction[XX];
        crntProbePos_[YY] = crntProbePos_[YY] + probeStepLength_*direction[YY];
        crntProbePos_[ZZ] = crntProbePos_[ZZ] + probeStepLength_*direction[ZZ]; 

        // optimise in plane through simulated annealing:
        SimulatedAnnealingModule sam;
        sam.setObjFun(objFun);
        sam.setParams(params_);
        sam.setInitGuess(initState);
        sam.optimise();

        // refine with Nelder-Mead optimisation:
        NelderMeadModule nmm;
        nmm.setObjFun(objFun);
        nmm.setParams(params_);
        nmm.setInitGuess(sam.getOptimPoint().first);
        nmm.optimise();
 
        // current position becomes best position in plane: 
        crntProbePos_ = optimToConfig(nmm.getOptimPoint().first);
               
        // increment probe step counter:
        numProbeSteps++;      

        // add result to path container: 
        path_.push_back(crntProbePos_);
        radii_.push_back(nmm.getOptimPoint().second);     

        // check termination conditions:
        if( numProbeSteps >= maxProbeSteps_ )
        {
            break;
        }
        if( nmm.getOptimPoint().second > maxProbeRadius_ )
        {
            break;
        }
    }

    // change radius of ultimate point to match the desired cutoff exactly:
    radii_.back() = maxProbeRadius_;
}


/*!
 * Converts between the two-dimensional optimisation space representation to 
 * the three-dimensional configuration space representation. A point in 
 * optimisation space is represented by its position on terms of the in-plane
 * basis spanning vectors orthVecU_ and orthVecW_, which are both orthogonal
 * to the channel direction vector.
 */
gmx::RVec
InplaneOptimisedProbePathFinder::optimToConfig(std::vector<real> optimSpacePos)
{
    // get configuration space position via orthogonal vectors:
    gmx::RVec configSpacePos;
    configSpacePos[XX] = crntProbePos_[XX] + optimSpacePos[0]*orthVecU_[XX]
                                           + optimSpacePos[1]*orthVecW_[XX];
    configSpacePos[YY] = crntProbePos_[YY] + optimSpacePos[0]*orthVecU_[YY]
                                           + optimSpacePos[1]*orthVecW_[YY];
    configSpacePos[ZZ] = crntProbePos_[ZZ] + optimSpacePos[0]*orthVecU_[ZZ] 
                                           + optimSpacePos[1]*orthVecW_[ZZ];
    
    // return configuration space position:
    return(configSpacePos);
}

