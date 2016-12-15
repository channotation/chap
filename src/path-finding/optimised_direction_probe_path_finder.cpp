#include <iostream>
#include <cmath>

#include "trajectoryAnalysis/simulated_annealing_module.hpp"
#include "path-finding/optimised_direction_probe_path_finder.hpp"

/*
 * Constructor.
 */
OptimisedDirectionProbePathFinder::OptimisedDirectionProbePathFinder(
        real probeStepLength,
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
        real saCoolingFactor,
        real saStepLengthFactor,
        bool saUseAdaptiveCandGen)
    : AbstractProbePathFinder(probeStepLength, probeRadius, maxFreeDist, 
                              maxProbeSteps, initProbePos, vdwRadii, nbSearch,
                              saRandomSeed, saMaxCoolingIter, saNumCostSamples,
                              saXi, saConvRelTol, saInitTemp, saCoolingFactor,
                              saStepLengthFactor, saUseAdaptiveCandGen)
{

    // initialise inverse rotation matrix as identity matrix:
    clear_mat(inverseRotationMatrix_);
}


/*
 *
 */
void
OptimisedDirectionProbePathFinder::findPath()
{
    std::cout<<"Yo!"<<std::endl;

    

}


/*
 *
 */
void
OptimisedDirectionProbePathFinder::optimiseInitialPos()
{

}


/*
 *
 */
void
OptimisedDirectionProbePathFinder::advanceAndOptimise(bool forward)
{
    gmx::RVec direction(0, 0, 1);
    if( !forward )
    {
        direction[0] = -direction[0];
        direction[1] = -direction[1];
        direction[2] = -direction[2];
    }

    // initial state in optimisation space is always the null vector:
    int stateDim = 2;
    real initState[stateDim] = {0.0, 0.0};

    // cost function is minimal free distance function:
    costFunction cf;
    cf = std::bind(&OptimisedDirectionProbePathFinder::findMinimalFreeDistance,
                   this, std::placeholders::_1);

    // advance probe position:
    int numProbeSteps = 0;
    while( true )
    {

        // create new simulated annealing module:                               
        SimulatedAnnealingModule sam(stateDim, saRandomSeed_, saMaxCoolingIter_,
                                     saNumCostSamples_, saXi_, saConvRelTol_,   
                                     saInitTemp_, saCoolingFactor_,                
                                     saStepLengthFactor_, initState, cf,           
                                     saUseAdaptiveCandGen_); 

        // optimise direction:
        eSimAnTerm status = sam.anneal();

        // opdate current and previous probe position:
        prevProbePos_ = crntProbePos_;
        crntProbePos_ = optimToConfig(sam.getBestState());

        // add result to path container:
        path_.push_back(crntProbePos_);
        radii_.push_back(sam.getBestCost());

        // update inverse rotation matrix:
        updateInverseRotationMatrix();

        // incremet probe step counter:
        numProbeSteps++;

        // termination conditions:
        if( sam.getBestCost() > maxProbeRadius_ )
        {
            break;
        }
        if( numProbeSteps >= maxProbeSteps_ )
        {
            break;
        }
    }
}


/*
 * Updates inverse rotation matrix used in conversion between optimisation and 
 * configuration space.
 */
void
OptimisedDirectionProbePathFinder::updateInverseRotationMatrix()
{
    // internal tolerance for floating point comparison:
    real zeroTol = 1e-7;

    // reset inverse rotation matrix to zero:
    clear_mat(inverseRotationMatrix_);

    // get z-basis vector in rotated and standard/global system:
    gmx::RVec stdBasisZ(0.0, 0.0, 1.0);
    gmx::RVec rotBasisZ;
    rvec_sub(crntProbePos_, prevProbePos_, rotBasisZ);
    unitv(rotBasisZ, rotBasisZ);

    // calculate cosine of rotation angle: 
    real cosRotAngle = iprod(stdBasisZ, rotBasisZ);

    // check if edge cases need to be handled:
    if( std::abs(cosRotAngle - 1.0) < zeroTol )
    {
        // in parallel case inverse rotation matrix is identity matrix:
        inverseRotationMatrix_[XX][XX] = 1.0;
        inverseRotationMatrix_[YY][YY] = 1.0;
        inverseRotationMatrix_[ZZ][ZZ] = 1.0;
    }
    else if (std::abs(cosRotAngle + 1.0) < zeroTol )
    {
        // in anti-parallel case rotate by pi around (global) x-axis:
        inverseRotationMatrix_[XX][XX] =  1.0;
        inverseRotationMatrix_[XX][XX] = -1.0;
        inverseRotationMatrix_[XX][XX] = -1.0; 
    }
    else
    {
        // calculate rotation axis vector:
        gmx::RVec rotAxisVec;
        cprod(stdBasisZ, rotBasisZ, rotAxisVec);

        // calculate cross product matrix:
        matrix rotationMatrix;
        rotationMatrix[XX][YY] = -rotAxisVec[ZZ];
        rotationMatrix[XX][ZZ] =  rotAxisVec[YY];
        rotationMatrix[YY][XX] =  rotAxisVec[ZZ];
        rotationMatrix[YY][ZZ] = -rotAxisVec[XX];
        rotationMatrix[ZZ][XX] = -rotAxisVec[YY];
        rotationMatrix[ZZ][YY] =  rotAxisVec[XX];

        // calculate and scale square cross product matrix:
        matrix sqCrossProdMat;
        mmul(rotationMatrix, rotationMatrix, sqCrossProdMat);
        msmul(sqCrossProdMat, 1.0/(1.0 + cosRotAngle), sqCrossProdMat);

        // add this to inverse rotation matrix:
        m_add(rotationMatrix, sqCrossProdMat, rotationMatrix);

        // add diagonal entries to inverse roation matrix:
        rotationMatrix[XX][XX] += 1.0;
        rotationMatrix[YY][YY] += 1.0;
        rotationMatrix[ZZ][ZZ] += 1.0;

        // transpose this:
        transpose(rotationMatrix, inverseRotationMatrix_);
    }    
}


/*
 * Performs optimisation space to configuration space conversion.
 *
 * The input array optimSpacePos contains the angles theta and phi in the local
 * (rotated) coordinate system that describe the probe direction vector:
 *
 * optimSpacePos[0] = theta
 * optimSpacePos[1] = phi
 *
 * This is first converted to a cartesian representation in the rotated system,
 * which is then mapped back to the global cartesian frame. The direction 
 * vector is then added to the position vector of the probe in the previous 
 * step to obtain the new configuration space position of the probe.
 */
gmx::RVec
OptimisedDirectionProbePathFinder::optimToConfig(real *optimSpacePos)
{
    // declare result vector:
    gmx::RVec configSpacePos;

    // calculate cartesian form of direction vector in rotated system:
    gmx::RVec rotatedDirVec(std::sin(optimSpacePos[0]) * std::cos(optimSpacePos[1]),
                            std::sin(optimSpacePos[0]) * std::sin(optimSpacePos[1]),
                            std::cos(optimSpacePos[0]));
    svmul(probeStepLength_, rotatedDirVec, rotatedDirVec);

    // rotate back to global system:
    mvmul(inverseRotationMatrix_, rotatedDirVec, configSpacePos);

    // add this to previous probe position:
    rvec_add(configSpacePos, crntProbePos_, configSpacePos);

    // return configuration space position:
    return(configSpacePos);
}

