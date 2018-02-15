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

#include <iostream>
#include <cmath>

#include "optim/simulated_annealing_module.hpp"
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

    // opimise initial
    
    // forward

    // array inversion

    // backward
}


/*
 *
 */
void
OptimisedDirectionProbePathFinder::optimiseInitialPos()
{

}


/*!
 * Advances the probe position in an optimised direction until the probe radius
 * exceeds a specified limit (or a maximum number of probe steps has been 
 * exceeded).
 */
void
OptimisedDirectionProbePathFinder::advanceAndOptimise(gmx::RVec initDirection)
{

}


/*!
 * Updates inverse rotation matrix used in conversion between optimisation and 
 * configuration space.
 *
 * The exact procedure is based on Rodrigues' rotation formula, from which it
 * can be derived that the rotation matrix mapping a unit vector a onto a unit
 * vector b is given by
 *
 *     R = I + K + K^2 * 1/(1+cos(alpha))
 *
 * Here I is the identity matrix, K is the cross-product matrix prescribing the
 * direction of the rotation axis in space, and alpha the angle by which 
 * vector a is rotated around that axis to be mapped onto b.
 */
void
OptimisedDirectionProbePathFinder::updateInverseRotationMatrix(gmx::RVec direction)
{
    /*

    // internal tolerance for floating point comparison:
    real zeroTol = 1e-7;

    // reset inverse rotation matrix to zero:
    clear_mat(inverseRotationMatrix_);

    // get z-basis vector in rotated and standard/global system:
    gmx::RVec stdBasisZ(0.0, 0.0, 1.0);
    gmx::RVec rotBasisZ(direction);
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

    */
}


/*!
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
OptimisedDirectionProbePathFinder::optimToConfig(std::vector<real> optimSpacePos)
{
    /*

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
    */
}

#endif

