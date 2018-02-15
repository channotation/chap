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

#include <gromacs/math/vec.h>

#include "optim/simulated_annealing_module.hpp"
#include "optim/nelder_mead_module.hpp"

#include "path-finding/inplane_optimised_probe_path_finder.hpp"



/*
 *
 */
InplaneOptimisedProbePathFinder::InplaneOptimisedProbePathFinder(
        std::map<std::string, real> params,
        gmx::RVec initProbePos,
        gmx::RVec chanDirVec,
//        gmx::AnalysisNeighborhoodSearch *nbSearch,
        t_pbc pbc,
        gmx::AnalysisNeighborhoodPositions porePos,
        std::vector<real> vdwRadii)
    : AbstractProbePathFinder(params, initProbePos, pbc, porePos, vdwRadii)
    , chanDirVec_(chanDirVec)
    , orthVecU_(0.0, 0.0, 0.0)
    , orthVecW_(0.0, 0.0, 0.0)
{
    // TODO: all these vector operations should be written as testable function!

    // normalise channel direction vector:
    unitv(chanDirVec_, chanDirVec_);

    // tolerance threshold for norm of vector (which should be unit vectors):
    real nonZeroTol = 1e-1;

    if( norm(chanDirVec_) < nonZeroTol )
    {
        std::cout<<"ERROR: channel direction vector is (almost) zero!"<<std::endl;
    }

    // generate first orthogonal vector:
    orthVecU_ = gmx::RVec(-chanDirVec_[1], chanDirVec_[0], 0.0);

    // make sure it is non null vector:
    if( norm(orthVecU_) < nonZeroTol )
    {
        // try different permutation:
        orthVecU_ = gmx::RVec(-chanDirVec_[2], 0.0, chanDirVec_[0]);

        // make sure it is not null vector:
        if( norm(orthVecU_) < nonZeroTol )
        {
            // try different permutation:
            orthVecU_ = gmx::RVec(0.0, -chanDirVec_[2], chanDirVec_[1]);

            // make sure it is not null vector:
            if( norm(orthVecU_) < nonZeroTol )
            {
                std::cerr<<"ERROR: could not generate orthogonal vector!"<<std::endl;
            }
        }
    }

    // generate second orthogonal vector:
    cprod(chanDirVec_, orthVecU_, orthVecW_);

    // make sure vctors are mutually orthogonal:
    // TODO: proper error handling!
    real orthTol = 1e-7;
    if( iprod(chanDirVec_, orthVecU_) > orthTol )
    {
        std::cout<<"ERROR: basis vectors not orthogonal!"<<std::endl;
    }
    if( iprod(chanDirVec_, orthVecW_) > orthTol )
    {
        std::cout<<"ERROR: basis vectors not orthogonal!"<<std::endl;
    }
    if( iprod(orthVecU_, orthVecW_) > orthTol )
    {
        std::cout<<"ERROR: basis vectors not orthogonal!"<<std::endl;
    }
 
}


/*
 *
 */
void
InplaneOptimisedProbePathFinder::findPath()
{
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


/*
 *
 */
void
InplaneOptimisedProbePathFinder::optimiseInitialPos()
{
    // set current probe position to initial probe position: 
    crntProbePos_ = initProbePos_;
/*
    std::cout<<"crntProbePos = "<<crntProbePos_[0]<<"  "
                                <<crntProbePos_[1]<<"  "
                                <<crntProbePos_[2]<<"  "
             <<"initProbePos = "<<initProbePos_[0]<<"  "
                                <<initProbePos_[1]<<"  "
                                <<initProbePos_[2]<<"  "
             <<"orthVecU = "<<orthVecU_[0]<<"  "
                            <<orthVecU_[1]<<"  "
                            <<orthVecU_[2]<<"  "
             <<"orthVecW = "<<orthVecW_[0]<<"  "
                            <<orthVecW_[1]<<"  "
                            <<orthVecW_[2]<<"  "
              <<std::endl;
*/


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


/*
    std::cout<<"initProbePos = "<<initProbePos_[0]<<"  "
                                <<initProbePos_[1]<<"  "
                                <<initProbePos_[2]<<"  "
             <<"status = "<<status<<"  "
             <<"bestState = "<<sam.getBestState()[0]<<" "
                             <<sam.getBestState()[1]<<" "
             <<std::endl;
    std::cout<<"crntProbePos = "<<crntProbePos_[0]<<"  "
                                <<crntProbePos_[1]<<"  "
                                <<crntProbePos_[2]<<"  "
             <<std::endl; 
*/


    // add path support point and associated radius to container:
    path_.push_back(initProbePos_);
    radii_.push_back(nmm.getOptimPoint().second);

   
}


/*
 *
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
        direction[0] = -direction[0];
        direction[1] = -direction[1];
        direction[2] = -direction[2];
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
    /*
        std::cout<<"probeStepLength = "<<probeStepLength_<<std::endl;
        std::cout<<"direction = "<<direction[0]<<"  "
                                 <<direction[1]<<"  "
                                 <<direction[2]<<"  "
                 <<"crntProbePos = "<<crntProbePos_[0]<<"  "
                                    <<crntProbePos_[1]<<"  "
                                    <<crntProbePos_[2]<<"  "
                 <<std::endl;
*/
        // advance probe position to next plane:
        crntProbePos_[0] = crntProbePos_[0] + probeStepLength_*direction[0];
        crntProbePos_[1] = crntProbePos_[1] + probeStepLength_*direction[1];
        crntProbePos_[2] = crntProbePos_[2] + probeStepLength_*direction[2]; 

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
       
        // add result to path container: 
        path_.push_back(crntProbePos_);
        radii_.push_back(nmm.getOptimPoint().second);     
        
        // increment probe step counter:
        numProbeSteps++;
//        std::cout<<"probe step = "<<numProbeSteps<<std::endl;
//        std::cout<<"sa best cost = "<<sam.getOptimPoint().second<<std::endl;
//        std::cout<<"nm best cost = "<<nmm.getOptimPoint().second<<std::endl;
//        std::cout<<"crntProbePos = "<<crntProbePos_[XX]<<"  "
//                                    <<crntProbePos_[YY]<<"  "
//                                    <<crntProbePos_[ZZ]<<"  "<<std::endl;

       
        if( numProbeSteps >= maxProbeSteps_ )
        {
            break;
        }
        if( nmm.getOptimPoint().second > maxProbeRadius_ )
        {
            break;
        }
    }
}


/*
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
    configSpacePos[0] = crntProbePos_[0] + optimSpacePos[0]*orthVecU_[0]
                                         + optimSpacePos[1]*orthVecW_[0];
    configSpacePos[1] = crntProbePos_[1] + optimSpacePos[0]*orthVecU_[1]
                                         + optimSpacePos[1]*orthVecW_[1];
    configSpacePos[2] = crntProbePos_[2] + optimSpacePos[0]*orthVecU_[2] 
                                         + optimSpacePos[1]*orthVecW_[2];
    
    // return configuration space position:
    return(configSpacePos);
}

