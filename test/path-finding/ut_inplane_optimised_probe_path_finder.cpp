#include <algorithm>
#include <functional>
#include <fstream>

#include <gtest/gtest.h>

#include "path-finding/inplane_optimised_probe_path_finder.hpp"
#include "path-finding/mock_pore_maker.hpp"


/*
 *
 */
class InplaneOptimisedProbePathFinderTest : public ::testing::Test
{

    public:

        // constructor:
        InplaneOptimisedProbePathFinderTest()
        {
                // path finder parameters:
                params_["pfProbeRadius"] = 0.0;
                params_["pfProbeStepLength"] = 0.01;
                params_["pfProbeMaxRadius"] = 1.0;
                params_["pfProbeMaxSteps"] = 1000;

                // simulated annealing parameters:
                params_["saUseAdaptiveCandidateGeneration"] = 0;
                params_["saRandomSeed"] = 15011991;
                params_["saMaxCoolingIter"] = 25;
                params_["saNumCostSamples"] = 10;
                params_["saXi"] = 3.0;
                params_["saConvRelTol"] = 1e-15;
                params_["saCoolingFactor"] = 0.98;
                params_["saInitTemp"] = 0.1;
                params_["saStepLengthFactor"] = 0.001;

                // Nelder-Mead parameters:
                params_["nmMaxIter"] = 25;
                params_["nmInitShift"] = 0.1;

                // set periodic boundary condition struct:
                // NOTE: box is chosen so that periodicity does not matter
                matrix boxMat;
                boxMat[XX][XX] = 100.0*params_["pfMaxProbeRadius"];
                boxMat[XX][YY] = 0.0;
                boxMat[XX][ZZ] = 0.0;
                boxMat[YY][XX] = 0.0;
                boxMat[YY][YY] = 100.0*params_["pfMaxProbeRadius"];
                boxMat[YY][ZZ] = 0.0;
                boxMat[ZZ][XX] = 0.0;
                boxMat[ZZ][YY] = 0.0;
                boxMat[ZZ][ZZ] = 100.0*params_["pfMaxProbeRadius"];
                set_pbc(&pbc_, -1, boxMat);
        };

        // standard parameters for tests:
        std::map<std::string, real> params_;

        // pbc struct:
        t_pbc pbc_;

        // comparison functions used in STL counting:
        static bool isGreaterThan(real arg, real lim)
        {
            return (arg > lim);
        };
        static bool isLesserThan(real arg, real lim)
        {
            return (arg < lim);
        };
};


/*
 *
 */
TEST_F(InplaneOptimisedProbePathFinderTest, InplaneOptimisedProbePathFinderZDirTest)
{
    // set parameters to defaults:
    std::map<std::string, real> params = params_;

    // define pore parameters:
    real poreLength = 2.0;
    real poreCentreRadius = 0.5;
    real poreVdwRadius = 0.2;
    gmx::RVec poreCentre(0.0, 0.0, 0.0);
    int poreDir = ZZ;

    // calculate true pore radius:
    real poreMinFreeRadius = poreCentreRadius - poreVdwRadius;
    real poreMaxFreeRadius = std::sqrt(std::pow(poreVdwRadius/4.0, 2.0) + 
                             std::pow(poreCentreRadius, 2.0)) - poreVdwRadius;

    // create pore pointing in the z-direction:
    MockPoreMaker mpm;
    std::vector<gmx::RVec> particleCentres = mpm.makePore(poreLength,
                                                          poreCentreRadius,
                                                          poreVdwRadius,
                                                          poreCentre,
                                                          poreDir);    
    std::vector<real> vdwRadii;
    vdwRadii.insert(vdwRadii.begin(), particleCentres.size(), poreVdwRadius);


    std::fstream file;
    file.open("mockpore.dat", std::ios::out);

    for(int i = 0; i < particleCentres.size(); i++)
    {
//        std::cout<<i<<std::endl;
        file<<particleCentres[i][XX]<<"\t"
            <<particleCentres[i][YY]<<"\t"
            <<particleCentres[i][ZZ]<<"\t"
            <<std::endl;
    }

    file.close();



    // prepare an analysis neighborhood:
    gmx::AnalysisNeighborhood nbh;
    nbh.setCutoff(0.0);
    nbh.setXYMode(false);
    nbh.setMode(gmx::AnalysisNeighborhood::eSearchMode_Grid);


    // prepare neighborhood search:
    gmx::AnalysisNeighborhoodPositions nbhPos(particleCentres);
    std::cout<<"nparticleCentres = "<<particleCentres.size()<<std::endl;
    std::cout<<"ePBC = "<<pbc_.ePBC<<std::endl;
    std::cout<<"ndim_ePBC = "<<pbc_.ndim_ePBC<<std::endl;
    std::cout<<"ePBCDX = "<<pbc_.ePBCDX<<std::endl;
    std::cout<<"dim = "<<pbc_.dim<<std::endl;
    std::cout<<"fbox_diag = "<<pbc_.fbox_diag[XX]<<"  "
                             <<pbc_.fbox_diag[YY]<<"  "
                             <<pbc_.fbox_diag[ZZ]<<std::endl;
    std::cout<<"hbox_diag = "<<pbc_.hbox_diag[XX]<<"  "
                             <<pbc_.hbox_diag[YY]<<"  "
                             <<pbc_.hbox_diag[ZZ]<<std::endl;
    std::cout<<"mhbox_diag = "<<pbc_.mhbox_diag[XX]<<"  "
                              <<pbc_.mhbox_diag[YY]<<"  "
                              <<pbc_.mhbox_diag[ZZ]<<std::endl;
    std::cout<<"max_cutoff2 = "<<pbc_.max_cutoff2<<std::endl;
    std::cout<<"ntric_vec = "<<pbc_.ntric_vec<<std::endl;
    std::cout<<"vdwRadii.size = "<<vdwRadii.size()<<std::endl;
    std::cout<<"particleCentres.size = "<<particleCentres.size()<<std::endl;
    std::cout<<"poreMinFreeRadius = "<<poreMinFreeRadius<<std::endl;

//    std::cout<<"ePBCDX = "<<ePBCDX<<std::endl;
//    std::cout<<"ePBCDX = "<<ePBCDX<<std::endl;
    
    gmx::AnalysisNeighborhoodSearch nbSearch = nbh.initSearch(&pbc_,
                                                              nbhPos);

    // prepare path finder:
    gmx::RVec initProbePos(poreCentre[XX] + 0.0*poreCentreRadius, 
                           poreCentre[YY] - 0.0*poreCentreRadius, 
                           poreCentre[ZZ] + 0.0*poreCentreRadius);
    gmx::RVec chanDirVec(0.0, 0.0, 1.0);

    // create path finder:
    InplaneOptimisedProbePathFinder pfm(params,
                                        initProbePos,
                                        chanDirVec,
//                                        &nbSearch,
                                        pbc_,
                                        nbhPos,
                                        vdwRadii);

    // find and extract path and path points:
    pfm.findPath();
//    MolecularPath molPath = pfm.getMolecularPath();
    std::vector<real> radii = pfm.pathRadii();
    std::vector<gmx::RVec> points = pfm.pathPoints();
    
    file.open("foundpore.dat", std::ios::out);

    for(int i = 0; i < points.size(); i++)
    {
//        std::cout<<i<<std::endl;
        file<<points[i][XX]<<"\t"
            <<points[i][YY]<<"\t"
            <<points[i][ZZ]<<"\t"
            <<radii[i]
            <<std::endl;
    }

    file.close();

    // check that no points are negative:
    std::function<bool(real)> lt;
    lt = std::bind(isLesserThan, std::placeholders::_1, 0.0);
    int nNegative = std::count_if(radii.begin(), radii.end(), lt);
    ASSERT_GE(0, nNegative);

    // check that no more than two points exceed the termination radius:
    std::function<bool(real)> gt;
    gt = std::bind(isGreaterThan, std::placeholders::_1, params["pfProbeMaxRadius"]);
    int nGreaterLimit = std::count_if(radii.begin(), radii.end(), gt);
    ASSERT_GE(2, nGreaterLimit);

    // extract the pore internal points:
    std::vector<real> internalRadii;
    std::vector<gmx::RVec> internalPoints;
    for(unsigned int i = 0; i < radii.size(); i++)
    {
        if( points[i][poreDir] >= poreCentre[poreDir] - 0.5*poreLength && 
            points[i][poreDir] <= poreCentre[poreDir] + 0.5*poreLength )
        {
            internalPoints.push_back(points[i]);
            internalRadii.push_back(radii[i]);
        }
    }

    // check that all internal radii have the minimal free distance:
    real minFreeRadTol = std::numeric_limits<real>::epsilon();
    for(unsigned int i = 0; i < internalRadii.size(); i++)
    {
        ASSERT_LE(poreMinFreeRadius - minFreeRadTol, internalRadii[i]);
    }

    // check that all internal radii have less than the maximal free distance:
    real maxFreeRadTol = std::numeric_limits<real>::epsilon();
    for(unsigned int i = 0; i < internalRadii.size(); i++)
    {
        std::cout<<"internalRad = "<<internalRadii[i]<<std::endl;
        std::cout<<"maxFreeRad = "<<poreMaxFreeRadius<<std::endl;

        ASSERT_GE(poreMaxFreeRadius + maxFreeRadTol, internalRadii[i]);
    }

    // check that all internal points lie on the centreline:
    real clDistTol = 10*std::numeric_limits<real>::epsilon();
    for(unsigned int i = 0; i < internalPoints.size(); i++)
    {
        // find local point on centreline:
        gmx::RVec centreLinePoint(poreCentre);
        centreLinePoint[poreDir] -= i*params["pfProbeStepLength"] - 0.5*poreLength;

        // assert vanishingb distance from centreline:
        real clDist = std::sqrt(distance2(internalPoints[i], centreLinePoint));
        ASSERT_NEAR(0.0, clDist, clDistTol);
    }


    file.open("foundporeinternal.dat", std::ios::out);

    for(int i = 0; i < internalPoints.size(); i++)
    {
//        std::cout<<i<<std::endl;
        file<<internalPoints[i][XX]<<"\t"
            <<internalPoints[i][YY]<<"\t"
            <<internalPoints[i][ZZ]<<"\t"
            <<internalRadii[i]
            <<std::endl;
    }

    file.close();





}



