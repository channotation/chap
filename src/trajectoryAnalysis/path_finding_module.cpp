#include <iostream>

#include "trajectoryAnalysis/path_finding_module.hpp"
#include "trajectoryAnalysis/simulated_annealing_module.hpp"


/*
 * Constructor.
 */
PathFindingModule::PathFindingModule(gmx::RVec initProbePos,
                                     gmx::RVec chanDirVec,
                                     gmx::AnalysisNeighborhoodSearch nbSearch,
                                     std::vector<real> vdwRadii)
    : stepLength_(0.5)
    , vdwRadii_(vdwRadii)
    , initProbePos_(initProbePos)
    , chanDirVec_(chanDirVec)
    , nbSearch_(nbSearch)
{
	// normalise direction vector:
    real norm = std::sqrt( chanDirVec_[0]*chanDirVec_[0] +
                           chanDirVec_[1]*chanDirVec_[1] +
                           chanDirVec_[2]*chanDirVec_[2]);
    chanDirVec_[0] = chanDirVec_[0]/norm;
    chanDirVec_[1] = chanDirVec_[1]/norm;
    chanDirVec_[2] = chanDirVec_[2]/norm;

    std::cout<<"initProbePos_ = "<<initProbePos_[0]<<"  "
                                 <<initProbePos_[1]<<"  "
                                 <<initProbePos_[2]<<std::endl;


    // set up neighbourhood search engine:
//    nbSearch_ = neighborhood.initSearch(pbc, poreSelection);
//    nbSearch_ = nbSearch;
}


/*
 * Destructor.
 */
PathFindingModule::~PathFindingModule()
{

}


/*
 *
 */
void
PathFindingModule::findPath()
{
    std::cout<<"initProbePos_ = "<<initProbePos_[0]<<"  "
                                 <<initProbePos_[1]<<"  "
                                 <<initProbePos_[2]<<std::endl;

    // positive channel direction:
    marchAndOptimise(initProbePos_, chanDirVec_);
    
    // negative nnel direction:
    // negChanDirVec
    // marchAndOptimise(initProbePos
}


/*
 *
 */
void
PathFindingModule::marchAndOptimise(gmx::RVec initPos,
                                    bool forward)
{   
    std::cout<<"initPos = "<<initPos[0]<<"  "
                           <<initPos[1]<<"  "
                           <<initPos[2]<<std::endl;

    // set previous position to initial point:
    crntProbePos_ = initPos;

    // set up direction vector for forward/backward marching:
    gmx::RVec direction(chanDirVec_);
    if( !forward )
    {
        direction[0] = -direction[0];
        direction[1] = -direction[1];
        direction[2] = -direction[2];
    }

    // create orthogonal vectors:
    orthVecU_[0] = 1;
    orthVecU_[1] = 0;
    orthVecU_[2] = 0;
    orthVecW_[0] = 0;
    orthVecW_[1] = 1;
    orthVecW_[2] = 0;


    // prepare simulated annealing module:
    int stateDim = 2;
    int randomSeed = 15011991;
    int maxCoolingIter = 1e3;
    int numCostSamples = 50;
    real xi = 3.0;
    real convRelTol = 1e-10;
    real initTemp = 10;
    real coolingFactor = 0.98;
    real simAnStepLengthFactor = 0.2;
    real initstate[stateDim] = {0.0, 0.0};
    bool useAdaptiveCandidateGeneration = false;
    costFunction cf = std::bind(&PathFindingModule::findMinimumDistance, this, std::placeholders::_1);
   
   
    // add initial position to path container:
    // TODO: remove this!
    path_.push_back(crntProbePos_);
    radii_.push_back(1.0);
 
    eSimAnTerm status = CONVERGENCE;
    int fudge = 0;
    while(true)
    {
        // advance probe position to next plane:
        crntProbePos_[0] = crntProbePos_[0] + stepLength_*direction[0];
        crntProbePos_[1] = crntProbePos_[1] + stepLength_*direction[1];
        crntProbePos_[2] = crntProbePos_[2] + stepLength_*direction[2]; 

        // create new simulated annealing module:
        SimulatedAnnealingModule sam(stateDim,
                                     randomSeed,
                                     maxCoolingIter,
                                     numCostSamples,
                                     xi,
                                     convRelTol,
                                     initTemp,
                                     coolingFactor,
                                     simAnStepLengthFactor,
                                     initstate,
                                     cf,
                                     useAdaptiveCandidateGeneration);

        // optimise in plane:
        status = sam.anneal();
          
        // current position becomes best position in plane: 
        crntProbePos_ = optimToConfig(sam.getBestState());

        std::cout<<"status = "<<status<<"  "
                 <<"bestCost = "<<sam.getBestCost()<<"  "
                 <<"bestState(0) = "<<sam.getBestStateAt(0)<<"  "
                 <<"bestState(1) = "<<sam.getBestStateAt(1)<<"  "
                 <<"x = "<<crntProbePos_[0]<<"  "
                 <<"y = "<<crntProbePos_[1]<<"  "
                 <<"z = "<<crntProbePos_[2]<<std::endl;

        // add result to path container: 
        path_.push_back(crntProbePos_);
        radii_.push_back(sam.getBestCost());
      
        
               


        fudge++;
        if(fudge>20)
        {
            break;
        }
    }
}


/*
 * Finds minimum distance to reference selection.
 */
real
PathFindingModule::findMinimumDistance(real *optimSpacePos)
{

    // internal variables:
    real pairDist;              // distance between probe and pore atom
    real poreAtomVdwRadius;     // van-der-Waals radius of pore atom
    // TODO: infinity
    real voidRadius=99;            // radius of maximal non-overlapping sphere

    // convert point in optimisation space to point in configuration space:
    gmx::AnalysisNeighborhoodPositions probePos(optimToConfig(optimSpacePos).as_vec());

    // prepare neighbourhood search:
    gmx::AnalysisNeighborhoodPairSearch nbPairSearch = nbSearch_.startPairSearch(probePos);
    gmx::AnalysisNeighborhoodPair pair;

    // loop over all pairs:
    while( nbPairSearch.findNextPair(&pair) )
    {
        // get pair distance:
        pairDist = std::sqrt(pair.distance2());

        // get vdW radius of reference atom:
        poreAtomVdwRadius = vdwRadii_.at(pair.refIndex());

        // update void radius if necessary:
        // TODO: factor in a probe radius!
        if( (pairDist - poreAtomVdwRadius) < voidRadius )
        {
//        std::cout<<"voidRadius = "<<voidRadius<<std::endl;
            voidRadius = pairDist - poreAtomVdwRadius;
        }
    }
 //       std::cout<<"voidRadius = "<<voidRadius<<std::endl;

    // return radius of maximal free sphere:
    return voidRadius;
}


/*
 * Converts point coordinates in optimisation space to point in configuration 
 * space.
 */
gmx::RVec
PathFindingModule::optimToConfig(real *optimSpacePos)
{
/*    std::cout<<"orthVecU = "<<orthVecU_[0]<<" "
                            <<orthVecU_[1]<<" "
                            <<orthVecU_[2]<<" "<<std::endl;
    std::cout<<"orthVecW = "<<orthVecW_[0]<<" "
                            <<orthVecW_[1]<<" "
                            <<orthVecW_[2]<<" "<<std::endl;

    std::cout<<"crntProbePos = "<<crntProbePos_[0]<<" "
                                <<crntProbePos_[1]<<" "
                                <<crntProbePos_[2]<<" "<<std::endl;

    std::cout<<"optimSpacePos = "<<optimSpacePos[0]<<" "
                                 <<optimSpacePos[1]<<" "<<std::endl;
*/
    gmx::RVec configSpacePos;
    configSpacePos[0] = crntProbePos_[0] + optimSpacePos[0]*orthVecU_[0]
                                         + optimSpacePos[1]*orthVecW_[0];
    configSpacePos[1] = crntProbePos_[1] + optimSpacePos[0]*orthVecU_[1]
                                         + optimSpacePos[1]*orthVecW_[1];
    configSpacePos[2] = crntProbePos_[2] + optimSpacePos[0]*orthVecU_[2] 
                                         + optimSpacePos[1]*orthVecW_[2];
/*
    std::cout<<"configProbePos = "<<configSpacePos[0]<<" "
                                  <<configSpacePos[1]<<" "
                                 <<configSpacePos[2]<<" "<<std::endl;
  */  
    return(configSpacePos);
}




















