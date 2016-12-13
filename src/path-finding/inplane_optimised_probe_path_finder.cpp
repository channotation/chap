#include <iostream>

#include "trajectoryAnalysis/simulated_annealing_module.hpp"
#include "path-finding/inplane_optimised_probe_path_finder.hpp"

/*
 * Constructor.
 */
InplaneOptimisedProbePathFinder::InplaneOptimisedProbePathFinder(
        real probeStepLength,
        gmx::RVec initProbePos,
        std::vector<real> vdwRadii,
        gmx::AnalysisNeighborhoodSearch nbSearch)
    : AbstractProbePathFinder(probeStepLength, initProbePos, vdwRadii, nbSearch)
    , chanDirVec_(0.0, 0.0, 1.0)
    , orthVecU_(1.0, 0.0, 0.0)
    , orthVecW_(0.0, 1.0, 0.0)
{
    // TODO: make sure vectors are unit vectors and orthogonal!    
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

    // initial state in optimisation space is always null vector:
    int stateDim = 2;
    real initState[stateDim] = {0.0, 0.0};
    
    // cost function is minimal free distance function:
    costFunction cf;
    cf = std::bind(&InplaneOptimisedProbePathFinder::findMinimalFreeDistance, 
                   this, std::placeholders::_1);
   
    // create new simulated annealing module:
    SimulatedAnnealingModule sam(stateDim, saRandomSeed_, saMaxCoolingIter_,
                                 saNumCostSamples_, saXi_, saConvRelTol_, 
                                 saInitTemp_, saCoolingFactor_, 
                                 saStepLengthFactor_, initState, cf,
                                 saUseAdaptiveCandGen_);
    
    // perform simulated annealing:
    eSimAnTerm status = sam.anneal();
    
    // set initial position to its optimal value:
    initProbePos_ = optimToConfig(sam.getBestState());


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

    // add path support point and associated radius to container:
    path_.push_back(initProbePos_);
    radii_.push_back(sam.getBestCost());
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
    int stateDim = 2;
    real initState[stateDim] = {0.0, 0.0};

    // cost function is minimal free distance function:
    costFunction cf;
    cf = std::bind(&InplaneOptimisedProbePathFinder::findMinimalFreeDistance, 
                   this, std::placeholders::_1);
     
 
    eSimAnTerm status;
    int fudge = 0;
    real maxRadius = 10.0f;
    while(true)
    {
        std::cout<<"probeStepLength = "<<probeStepLength_<<std::endl;
        std::cout<<"direction = "<<direction[0]<<"  "
                                 <<direction[1]<<"  "
                                 <<direction[2]<<"  "
                 <<"crntProbePos = "<<crntProbePos_[0]<<"  "
                                    <<crntProbePos_[1]<<"  "
                                    <<crntProbePos_[2]<<"  "
                 <<std::endl;

        // advance probe position to next plane:
        crntProbePos_[0] = crntProbePos_[0] + probeStepLength_*direction[0];
        crntProbePos_[1] = crntProbePos_[1] + probeStepLength_*direction[1];
        crntProbePos_[2] = crntProbePos_[2] + probeStepLength_*direction[2]; 

        // create new simulated annealing module:
        SimulatedAnnealingModule sam(stateDim, saRandomSeed_, saMaxCoolingIter_,
                                     saNumCostSamples_, saXi_, saConvRelTol_,
                                     saInitTemp_, saCoolingFactor_,
                                     saStepLengthFactor_, initState, cf,
                                     saUseAdaptiveCandGen_);

        // optimise in plane:
        status = sam.anneal();
          
        // current position becomes best position in plane: 
        crntProbePos_ = optimToConfig(sam.getBestState());
        
        // add result to path container: 
        path_.push_back(crntProbePos_);
        radii_.push_back(sam.getBestCost());     
              
       


        fudge++;
        if(fudge>15)
        {
            break;
        }
        if( sam.getBestCost() > maxRadius )
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
InplaneOptimisedProbePathFinder::optimToConfig(real *optimSpacePos)
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
