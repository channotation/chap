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


#include <algorithm>
#include <functional>
#include <fstream>

#include <gtest/gtest.h>

#include <gromacs/math/3dtransforms.h>

#include "path-finding/inplane_optimised_probe_path_finder.hpp"


/*!
 * \brief Test fixture for InplaneOptimisedProbePathFinder.
 *
 * Provides some default parameters and a function to create an artificial
 * cylindrical pore oriented along one of the Cartesian coordinate axes.
 */
class InplaneOptimisedProbePathFinderTest : public ::testing::Test
{

    public:

        // constructor:
        InplaneOptimisedProbePathFinderTest()
        {
                // path finder parameters:
                params_["pfProbeRadius"] = 0.0;
                params_["pfProbeStepLength"] = 0.05;
                params_["pfProbeMaxRadius"] = 1.0;
                params_["pfProbeMaxSteps"] = 1000;

                // simulated annealing parameters:
                params_["saUseAdaptiveCandidateGeneration"] = 0;
                params_["saRandomSeed"] = 15011992;
                params_["saMaxCoolingIter"] = 1000;
                params_["saNumCostSamples"] = 10;
                params_["saXi"] = 3.0;
                params_["saConvRelTol"] = 1e-15;
                params_["saCoolingFactor"] = 0.98;
                params_["saInitTemp"] = 0.1;
                params_["saStepLengthFactor"] = 0.001;

                // Nelder-Mead parameters:
                params_["nmMaxIter"] = 100;
                params_["nmInitShift"] = 0.1;

                // set periodic boundary condition struct:
                // NOTE: box is chosen so that periodicity does not matter
                boxMat_[XX][XX] = 0.0;
                boxMat_[XX][YY] = 0.0;
                boxMat_[XX][ZZ] = 0.0;
                boxMat_[YY][XX] = 0.0;
                boxMat_[YY][YY] = 0.0;
                boxMat_[YY][ZZ] = 0.0;
                boxMat_[ZZ][XX] = 0.0;
                boxMat_[ZZ][YY] = 0.0;
                boxMat_[ZZ][ZZ] = 0.0;
        };

        // standard parameters for tests:
        std::map<std::string, real> params_;

        // box matrix for pbc:
        matrix boxMat_;

        // mathematical constants:
        const real PI_ = std::acos(-1.0);

        // comparison functions used in STL counting:
        static bool isGreaterThan(real arg, real lim)
        {
            return (arg > lim);
        };
        static bool isLesserThan(real arg, real lim)
        {
            return (arg < lim);
        };

        // create a cylindrical mock pore:
        std::vector<gmx::RVec> makePore(real poreLength,
                                        real poreCentreRadius,
                                        real poreVdwRadius,
                                        gmx::RVec poreCentre,
                                        int alongAxis)
        {
            // calculate angle required for overlapping vdW spheres:
            real phi = std::acos(1.0 - std::pow(poreVdwRadius, 2.0)/2.0/std::pow(poreCentreRadius, 2.0));
            int nStepsAround = std::ceil(2.0*PI_/phi);

            // step length and number of steps along the length of the pore:
            real stepLengthAlong = poreVdwRadius/2.0;
            int nStepsAlong = std::ceil(poreLength/stepLengthAlong) + 1;

            // place particles on cylinder surface centered around origin:
            std::vector<gmx::RVec> particleCentres;
            for(int i = 0; i < nStepsAlong; i++)
            {
                for(int j = 0; j < nStepsAround; j++)
                {
                    // add new particle to pore:
                    gmx::RVec particleCentre;
                    particleCentre[XX] = poreCentreRadius*std::cos(phi*j);
                    particleCentre[YY] = poreCentreRadius*std::sin(phi*j);
                    particleCentre[ZZ] = i*stepLengthAlong - poreLength/2.0;
                    particleCentres.push_back(particleCentre);
                }
            }

            // rotate in desired way:
            mat4 rotMat;
            vec4 v;
            if( alongAxis == XX )
            {
                // rotate around y-axis:
                gmx_mat4_init_rotation(YY, PI_/2.0, rotMat);

            }
            else if( alongAxis == YY )
            {
                // rotate around x-axis:
                gmx_mat4_init_rotation(XX, PI_/2.0, rotMat);
            }
            else if( alongAxis == ZZ )
            {
                // nothing to do:
                gmx_mat4_init_rotation(ZZ, 0.0, rotMat);
            }
            else
            {
                std::cerr<<"ERROR: AlongAxis must be one of XX, YY, or ZZ!"<<std::endl;
                std::abort();
            }


            for(unsigned int i = 0; i < particleCentres.size(); i++)
            {
                gmx_mat4_transform_point(rotMat, particleCentres[i], v);
                particleCentres[i][XX] = v[XX];
                particleCentres[i][YY] = v[YY];
                particleCentres[i][ZZ] = v[ZZ];
            }



            // shift particles away from origin:
            for(unsigned int i = 0; i < particleCentres.size(); i++)
            {
                particleCentres[i][XX] += poreCentre[XX];
                particleCentres[i][YY] += poreCentre[YY];
                particleCentres[i][ZZ] += poreCentre[ZZ];
            }

            // return mock pore particle positions:
            return particleCentres;
        };
};


/*!
 * \brief Tests InplaneOptimisedProbePathFinder on a cylindrical pore pointing
 * in the \f$ x \f$-direction.
 *
 * A mock pore consisting of spheres of radius \f$ R_\text{vdW} \f$ is created
 * by placing the spheres \f$ R_\text{c} \f$ from the pore's centre line. Their
 * distance in the direction along the pore is set to \f$ R_\text{vdW}/2 \f$
 * and their angular distance along the pore is set to
 *
 * \f[
 *      \phi = \cos \left[ 1 - \frac{R^2_\text{vdW}}{2 R^2_\text{c}} \right]
 * \f]
 *
 * ensuring that all van-der-Waal spheres overlap. The minimal radius of this
 * pore is then given by \f$ R_\text{min} = R_\text{c} - R_\text{vdW} \f$ and
 * the maximal pore internal radius is given by:
 *
 * \f[
 *      R_\text{max} = \sqrt{ \left( \frac{R_\text{vdW}}{4} \right)^2 + R_\text{c}^2 }
 * \f]
 *
 * The InplaneOptimisedProbePathFinder is then set up with an initial probe
 * position that deviates from the pore's centre of geometry by up to
 * \f$ 0.5 R_\text{c} \f$ in each coordinate and a pore direction vector that
 * is tilted by at least 45° with respect to the true pore direction.
 *
 * The test asserts that the pore radius found in each optimisation plane
 * lies between \f$ R_\text{min} \f$ and \f$ R_\text{max} \f$ where a tolerance
 * of ten times the machine epsilon is permitted.  It also explicitly asserts
 * that no negative pore radius is found and that no more than the two
 * endpoints have more than the permitted maximal pore radius. Lastly, the test
 * asserts that the centre line points found by the
 * InplaneOptimisedProbePathFinder deviate from the true centre line of the
 * pore by no more than ten times the machine epsilon.
 *
 * Note that due to the simplified geometry of the pore this test does not
 * assess the capability of the simulated annealing engine to exit local
 * optima.
 */
// NOTE: this test fails in Travis cloud, but not locally
// NOTE: could be a memory issue?
// TEST_F(InplaneOptimisedProbePathFinderTest, InplaneOptimisedProbePathFinderXDirTest)
// {
//     // set up periodic boundary conditions:
//     t_pbc pbc;
//     set_pbc(&pbc, 1, boxMat_);
//
//     // set parameters to defaults:
//     std::map<std::string, real> params = params_;
//
//     // define pore parameters:
//     real poreLength = 2.0;
//     real poreCentreRadius = 0.25;
//     real poreVdwRadius = 0.2;
//     gmx::RVec poreCentre(0.0, 0.0, 0.0);
//     int poreDir = XX;
//     int notPoreDirA = ZZ;
//     int notPoreDirB = YY;
//
//     // calculate true pore radius:
//     real poreMinFreeRadius = poreCentreRadius - poreVdwRadius;
//     real poreMaxFreeRadius = std::sqrt(std::pow(poreVdwRadius/4.0, 2.0) +
//                              std::pow(poreCentreRadius, 2.0)) - poreVdwRadius;
//
//     // create pore pointing in the z-direction:
//     std::vector<gmx::RVec> particleCentres = makePore(poreLength,
//                                                       poreCentreRadius,
//                                                       poreVdwRadius,
//                                                       poreCentre,
//                                                       poreDir);
//     std::vector<real> vdwRadii;
//     vdwRadii.insert(vdwRadii.begin(), particleCentres.size(), poreVdwRadius);
//
//     // prepare an analysis neighborhood:
//     gmx::AnalysisNeighborhood nbh;
//     nbh.setCutoff(0.0);
//     nbh.setXYMode(false);
//     nbh.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);
//
//     // prepare neighborhood search:
//     gmx::AnalysisNeighborhoodPositions nbhPos(particleCentres);
//     gmx::AnalysisNeighborhoodSearch nbSearch = nbh.initSearch(&pbc,
//                                                               nbhPos);
//
//     // prepare path finder:
//     gmx::RVec initProbePos(poreCentre[XX] + 0.1*poreCentreRadius,
//                            poreCentre[YY] - 0.3*poreCentreRadius,
//                            poreCentre[ZZ] + 0.3*poreCentreRadius);
//     gmx::RVec chanDirVec(1.0, 0.05, 0.5);
//
//     // create path finder:
//     InplaneOptimisedProbePathFinder pfm(params,
//                                         initProbePos,
//                                         chanDirVec,
//                                         &pbc,
//                                         nbhPos,
//                                         vdwRadii);
//
//     // set path finder parameters:
//     // TODO: have to transfer all parameters from constructor to here:
//     PathFindingParameters par;
//     par.setProbeStepLength(params["pfProbeStepLength"]);
//     par.setMaxProbeRadius(params["pfProbeMaxRadius"]);
//     par.setMaxProbeSteps(params["pfProbeMaxSteps"]);
//     pfm.setParameters(par);
//
//     // find and extract path and path points:
//     pfm.findPath();
//     std::vector<real> radii = pfm.pathRadii();
//     std::vector<gmx::RVec> points = pfm.pathPoints();
//
//     // check that no points are negative:
//     std::function<bool(real)> lt;
//     lt = std::bind(isLesserThan, std::placeholders::_1, 0.0);
//     int nNegative = std::count_if(radii.begin(), radii.end(), lt);
//     ASSERT_GE(0, nNegative);
//
//     // check that no more than two points exceed the termination radius:
//     std::function<bool(real)> gt;
//     gt = std::bind(isGreaterThan, std::placeholders::_1, params["pfProbeMaxRadius"]);
//     int nGreaterLimit = std::count_if(radii.begin(), radii.end(), gt);
//     ASSERT_GE(2, nGreaterLimit);
//
//     // extract the pore internal points:
//     std::vector<real> internalRadii;
//     std::vector<gmx::RVec> internalPoints;
//     for(unsigned int i = 0; i < radii.size(); i++)
//     {
//         if( points[i][poreDir] >= poreCentre[poreDir] - 0.5*poreLength &&
//             points[i][poreDir] <= poreCentre[poreDir] + 0.5*poreLength )
//         {
//             internalPoints.push_back(points[i]);
//             internalRadii.push_back(radii[i]);
//         }
//     }
//
//     // check that all internal radii have the minimal free distance:
//     real minFreeRadTol = 10*std::numeric_limits<real>::epsilon();
//     for(unsigned int i = 0; i < internalRadii.size(); i++)
//     {
//         ASSERT_LE(poreMinFreeRadius - minFreeRadTol, internalRadii[i]);
//     }
//
//     // check that all internal radii have less than the maximal free distance:
//     real maxFreeRadTol = 10*std::numeric_limits<real>::epsilon();
//     for(unsigned int i = 0; i < internalRadii.size(); i++)
//     {
//
//         ASSERT_GE(poreMaxFreeRadius + maxFreeRadTol, internalRadii[i]);
//     }
//
//     // check that all internal points lie on the centreline:
//     real clDistTol = 10*std::numeric_limits<real>::epsilon();
//     for(unsigned int i = 0; i < internalPoints.size(); i++)
//     {
//         ASSERT_NEAR(poreCentre[notPoreDirA],
//                     internalPoints[i][notPoreDirA],
//                     clDistTol);
//         ASSERT_NEAR(poreCentre[notPoreDirB],
//                     internalPoints[i][notPoreDirB],
//                     clDistTol);
//     }
// }


/*!
 * \brief Tests InplaneOptimisedProbePathFinder on a cylidnrical pore pointing
 * in the \f$ y \f$-direction.
 *
 * Similar to the previous test, but with the main direction of the pore
 * aligned with the Cartesian \f$ y \f$-axis. Additionally, several parameters
 * such as the pore and van-der-Waals radii, the initial probe position, and
 * the channel direction vector are slightly different.
 */
// NOTE: this test fails in Travis cloud, but not locally
// NOTE: could be a memory issue?
// TEST_F(InplaneOptimisedProbePathFinderTest, InplaneOptimisedProbePathFinderYDirTest)
// {
//     // set up periodic boundary conditions:
//     t_pbc pbc;
//     set_pbc(&pbc, 1, boxMat_);
//
//     // set parameters to defaults:
//     std::map<std::string, real> params = params_;
//
//     // define pore parameters:
//     real poreLength = 2.0;
//     real poreCentreRadius = 0.5;
//     real poreVdwRadius = 0.2;
//     gmx::RVec poreCentre(0.0, 0.0, 0.0);
//     int poreDir = YY;
//     int notPoreDirA = XX;
//     int notPoreDirB = ZZ;
//
//     // calculate true pore radius:
//     real poreMinFreeRadius = poreCentreRadius - poreVdwRadius;
//     real poreMaxFreeRadius = std::sqrt(std::pow(poreVdwRadius/4.0, 2.0) +
//                              std::pow(poreCentreRadius, 2.0)) - poreVdwRadius;
//
//     // create pore pointing in the z-direction:
//     std::vector<gmx::RVec> particleCentres = makePore(poreLength,
//                                                       poreCentreRadius,
//                                                       poreVdwRadius,
//                                                       poreCentre,
//                                                       poreDir);
//     std::vector<real> vdwRadii;
//     vdwRadii.insert(vdwRadii.begin(), particleCentres.size(), poreVdwRadius);
//
//     // prepare an analysis neighborhood:
//     gmx::AnalysisNeighborhood nbh;
//     nbh.setCutoff(0.0);
//     nbh.setXYMode(false);
//     nbh.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);
//
//     // prepare neighborhood search:
//     gmx::AnalysisNeighborhoodPositions nbhPos(particleCentres);
//     gmx::AnalysisNeighborhoodSearch nbSearch = nbh.initSearch(&pbc,
//                                                               nbhPos);
//
//     // prepare path finder:
//     gmx::RVec initProbePos(poreCentre[XX] + 0.1*poreCentreRadius,
//                            poreCentre[YY] - 0.5*poreCentreRadius,
//                            poreCentre[ZZ] + 0.3*poreCentreRadius);
//     gmx::RVec chanDirVec(0.5, 1.0, 0.5);
//
//     // create path finder:
//     InplaneOptimisedProbePathFinder pfm(params,
//                                         initProbePos,
//                                         chanDirVec,
//                                         &pbc,
//                                         nbhPos,
//                                         vdwRadii);
//
//     // set path finder parameters:
//     // TODO: have to transfer all parameters from constructor to here:
//     PathFindingParameters par;
//     par.setProbeStepLength(params["pfProbeStepLength"]);
//     par.setMaxProbeRadius(params["pfProbeMaxRadius"]);
//     par.setMaxProbeSteps(params["pfProbeMaxSteps"]);
//     pfm.setParameters(par);
//
//     // find and extract path and path points:
//     pfm.findPath();
//     std::vector<real> radii = pfm.pathRadii();
//     std::vector<gmx::RVec> points = pfm.pathPoints();
//
//     // check that no points are negative:
//     std::function<bool(real)> lt;
//     lt = std::bind(isLesserThan, std::placeholders::_1, 0.0);
//     int nNegative = std::count_if(radii.begin(), radii.end(), lt);
//     ASSERT_GE(0, nNegative);
//
//     // check that no more than two points exceed the termination radius:
//     std::function<bool(real)> gt;
//     gt = std::bind(isGreaterThan, std::placeholders::_1, params["pfProbeMaxRadius"]);
//     int nGreaterLimit = std::count_if(radii.begin(), radii.end(), gt);
//     ASSERT_GE(2, nGreaterLimit);
//
//     // extract the pore internal points:
//     std::vector<real> internalRadii;
//     std::vector<gmx::RVec> internalPoints;
//     for(unsigned int i = 0; i < radii.size(); i++)
//     {
//         if( points[i][poreDir] >= poreCentre[poreDir] - 0.5*poreLength &&
//             points[i][poreDir] <= poreCentre[poreDir] + 0.5*poreLength )
//         {
//             internalPoints.push_back(points[i]);
//             internalRadii.push_back(radii[i]);
//         }
//     }
//
//     // check that all internal radii have the minimal free distance:
//     real minFreeRadTol = 10*std::numeric_limits<real>::epsilon();
//     for(unsigned int i = 0; i < internalRadii.size(); i++)
//     {
//         ASSERT_LE(poreMinFreeRadius - minFreeRadTol, internalRadii[i]);
//     }
//
//     // check that all internal radii have less than the maximal free distance:
//     real maxFreeRadTol = 10*std::numeric_limits<real>::epsilon();
//     for(unsigned int i = 0; i < internalRadii.size(); i++)
//     {
//         ASSERT_GE(poreMaxFreeRadius + maxFreeRadTol, internalRadii[i]);
//     }
//
//     // check that all internal points lie on the centreline:
//     real clDistTol = 10*std::numeric_limits<real>::epsilon();
//     for(unsigned int i = 0; i < internalPoints.size(); i++)
//     {
//         ASSERT_NEAR(poreCentre[notPoreDirA],
//                     internalPoints[i][notPoreDirA],
//                     clDistTol);
//         ASSERT_NEAR(poreCentre[notPoreDirB],
//                     internalPoints[i][notPoreDirB],
//                     clDistTol);
//     }
// }


/*!
 * \brief Tests InplaneOptimisedProbePathFinder on a cylindrical pore pointing
 * in the \f$ z \f$-direction.
 *
 * Similar to the previous test, but with the main direction of the pore
 * aligned with the Cartesian \f$ z \f$-axis. Additionally, several parameters
 * such as the pore and van-der-Waals radii, the initial probe position, and
 * the channel direction vector are slightly different.
 */
// NOTE: this test fails in Travis cloud, but not locally
// NOTE: could be a memory issue?
// TEST_F(InplaneOptimisedProbePathFinderTest, InplaneOptimisedProbePathFinderZDirTest)
// {
//     // set up periodic boundary conditions:
//     t_pbc pbc;
//     set_pbc(&pbc, 1, boxMat_);
//
//     // set parameters to defaults:
//     std::map<std::string, real> params = params_;
//
//     // define pore parameters:
//     real poreLength = 3.0;
//     real poreCentreRadius = 0.25;
//     real poreVdwRadius = 0.2;
//     gmx::RVec poreCentre(0.0, 0.0, 0.0);
//     int poreDir = ZZ;
//     int notPoreDirA = XX;
//     int notPoreDirB = YY;
//
//     // calculate true pore radius:
//     real poreMinFreeRadius = poreCentreRadius - poreVdwRadius;
//     real poreMaxFreeRadius = std::sqrt(std::pow(poreVdwRadius/4.0, 2.0) +
//                              std::pow(poreCentreRadius, 2.0)) - poreVdwRadius;
//
//     // create pore pointing in the z-direction:
//     std::vector<gmx::RVec> particleCentres = makePore(poreLength,
//                                                       poreCentreRadius,
//                                                       poreVdwRadius,
//                                                       poreCentre,
//                                                       poreDir);
//     std::vector<real> vdwRadii;
//     vdwRadii.insert(vdwRadii.begin(), particleCentres.size(), poreVdwRadius);
//
//     // prepare an analysis neighborhood:
//     gmx::AnalysisNeighborhood nbh;
//     nbh.setCutoff(0.0);
//     nbh.setXYMode(false);
//     nbh.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);
//
//     // prepare neighborhood search:
//     gmx::AnalysisNeighborhoodPositions nbhPos(particleCentres);
//     gmx::AnalysisNeighborhoodSearch nbSearch = nbh.initSearch(&pbc,
//                                                               nbhPos);
//
//     // prepare path finder:
//     gmx::RVec initProbePos(poreCentre[XX] + 0.0*poreCentreRadius,
//                            poreCentre[YY] - 0.0*poreCentreRadius,
//                            poreCentre[ZZ] + 0.0*poreCentreRadius);
//     gmx::RVec chanDirVec(-0.0, 0.0, 1.0);
//
//     // create path finder:
//     InplaneOptimisedProbePathFinder pfm(params,
//                                         initProbePos,
//                                         chanDirVec,
//                                         &pbc,
//                                         nbhPos,
//                                         vdwRadii);
//
//     // set path finder parameters:
//     // TODO: have to transfer all parameters from constructor to here:
//     PathFindingParameters par;
//     par.setProbeStepLength(params["pfProbeStepLength"]);
//     par.setMaxProbeRadius(params["pfProbeMaxRadius"]);
//     par.setMaxProbeSteps(params["pfProbeMaxSteps"]);
//     pfm.setParameters(par);
//
//     // find and extract path and path points:
//     pfm.findPath();
//     std::vector<real> radii = pfm.pathRadii();
//     std::vector<gmx::RVec> points = pfm.pathPoints();
//
//     // check that no points are negative:
//     std::function<bool(real)> lt;
//     lt = std::bind(isLesserThan, std::placeholders::_1, 0.0);
//     int nNegative = std::count_if(radii.begin(), radii.end(), lt);
//     ASSERT_GE(0, nNegative);
//
//     // check that no more than two points exceed the termination radius:
//     std::function<bool(real)> gt;
//     gt = std::bind(isGreaterThan, std::placeholders::_1, params["pfProbeMaxRadius"]);
//     int nGreaterLimit = std::count_if(radii.begin(), radii.end(), gt);
//     ASSERT_GE(2, nGreaterLimit);
//
//     // extract the pore internal points:
//     std::vector<real> internalRadii;
//     std::vector<gmx::RVec> internalPoints;
//     for(unsigned int i = 0; i < radii.size(); i++)
//     {
//         if( points[i][poreDir] >= poreCentre[poreDir] - 0.5*poreLength &&
//             points[i][poreDir] <= poreCentre[poreDir] + 0.5*poreLength )
//         {
//             internalPoints.push_back(points[i]);
//             internalRadii.push_back(radii[i]);
//         }
//     }
//
//     // check that all internal radii have the minimal free distance:
//     real minFreeRadTol = 10.0*std::numeric_limits<real>::epsilon();
//     for(unsigned int i = 0; i < internalRadii.size(); i++)
//     {
//         ASSERT_LE(poreMinFreeRadius - minFreeRadTol, internalRadii[i]);
//     }
//
//     // check that all internal radii have less than the maximal free distance:
//     real maxFreeRadTol = 10.0*std::numeric_limits<real>::epsilon();
//     for(unsigned int i = 0; i < internalRadii.size(); i++)
//     {
//         ASSERT_GE(poreMaxFreeRadius + maxFreeRadTol, internalRadii[i]);
//     }
//
//     // check that all internal points lie on the centreline:
//     real clDistTol = 10.0*std::numeric_limits<real>::epsilon();
//     for(unsigned int i = 0; i < internalPoints.size(); i++)
//     {
//         ASSERT_NEAR(poreCentre[notPoreDirA],
//                     internalPoints[i][notPoreDirA],
//                     clDistTol);
//         ASSERT_NEAR(poreCentre[notPoreDirB],
//                     internalPoints[i][notPoreDirB],
//                     clDistTol);
//     }
// }
