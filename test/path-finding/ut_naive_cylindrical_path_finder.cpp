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


#include <gtest/gtest.h>

#include "path-finding/naive_cylindrical_path_finder.hpp"


/*!
 * \brief Test fixture for NaiveCylindricalPathFinder.
 *
 * Provides no functionality except that inherited from Google test.
 */
class NaiveCylindricalPathFinderTest : public ::testing::Test
{

    public:


};


/*!
 * \brief Tests the correct construction of a cylinder aligned with the 
 * \f$ x \f$-axis.
 *
 * The test asserts that both the centre line points and path radii are correct
 * to within one machine epsilon.
 */
TEST_F(NaiveCylindricalPathFinderTest, NaiveCylindricalPathFinderXDirTest)
{
    // set parameters:
    std::map<std::string, real> params;
    params["pfCylNumSteps"] = 10;
    params["pfCylStepLength"] = 0.01;
    params["pfCylRad"] = 1;

    // set direction vector and centre point:
    gmx::RVec dirVec(1.0, 0.0, 0.0);
    gmx::RVec centrePoint(0.0, 0.0, 0.0);

    // create path finding module and create path:
    NaiveCylindricalPathFinder pfm(params, centrePoint, dirVec);
    pfm.findPath();

    // find path:
    MolecularPath cylPath = pfm.getMolecularPath();

    // get path radii and support points:
    std::vector<real> radii = cylPath.pathRadii();
    std::vector<gmx::RVec> points = cylPath.pathPoints();

    // assert correct number of points:
    ASSERT_EQ(2*params["pfCylNumSteps"] + 1, radii.size());

    // assert correct radii:
    for(unsigned int i = 0; i < radii.size(); i++)
    {
        ASSERT_NEAR(params["pfCylRad"],
                    radii[i],
                    std::numeric_limits<real>::epsilon());
    }

    // assert correct location of path points:
    for(unsigned int i = 0; i < radii.size(); i++)
    {   
        // x-coordinate:
        ASSERT_NEAR(centrePoint[XX] + (i - params["pfCylNumSteps"])*params["pfCylStepLength"],
                    points[i][XX],
                    std::numeric_limits<real>::epsilon());

        // y-coordinate:
        ASSERT_NEAR(centrePoint[YY],
                    points[i][YY],
                    std::numeric_limits<real>::epsilon());

        // z-coordinate:
        ASSERT_NEAR(centrePoint[ZZ],
                    points[i][ZZ],
                    std::numeric_limits<real>::epsilon());
    }
}


/*!
 * \brief Tests the correct construction of a cylinder aligned with the 
 * \f$ y \f$-axis.
 *
 * The test asserts that both the centre line points and path radii are correct
 * to within one machine epsilon.
 */
TEST_F(NaiveCylindricalPathFinderTest, NaiveCylindricalPathFinderYDirTest)
{
    // set parameters:
    std::map<std::string, real> params;
    params["pfCylNumSteps"] = 5;
    params["pfCylStepLength"] = 1.0;
    params["pfCylRad"] = 0.025;

    // set direction vector and centre point:
    gmx::RVec dirVec(0.0, 1.0, 0.0);
    gmx::RVec centrePoint(4.5, -0.125, std::sqrt(2.0));

    // create path finding module and create path:
    NaiveCylindricalPathFinder pfm(params, centrePoint, dirVec);
    pfm.findPath();

    // find path:
    MolecularPath cylPath = pfm.getMolecularPath();

    // get path radii and support points:
    std::vector<real> radii = cylPath.pathRadii();
    std::vector<gmx::RVec> points = cylPath.pathPoints();

    // assert correct number of points:
    ASSERT_EQ(2*params["pfCylNumSteps"] + 1, radii.size());

    // assert correct radii:
    for(unsigned int i = 0; i < radii.size(); i++)
    {
        ASSERT_NEAR(params["pfCylRad"],
                    radii[i],
                    std::numeric_limits<real>::epsilon());
    }

    // assert correct location of path points:
    for(unsigned int i = 0; i < radii.size(); i++)
    {
        // x-coordinate:
        ASSERT_NEAR(centrePoint[XX],
                    points[i][XX],
                    std::numeric_limits<real>::epsilon());

        // y-coordinate:
        ASSERT_NEAR(centrePoint[YY] + (i - params["pfCylNumSteps"])*params["pfCylStepLength"],
                    points[i][YY],
                    std::numeric_limits<real>::epsilon());

        // z-coordinate:
        ASSERT_NEAR(centrePoint[ZZ],
                    points[i][ZZ],
                    std::numeric_limits<real>::epsilon());
    }
}


/*!
 * \brief Tests the correct construction of a cylinder aligned with the 
 * \f$ z \f$-axis.
 *
 * The test asserts that both the centre line points and path radii are correct
 * to within one machine epsilon.
 */
TEST_F(NaiveCylindricalPathFinderTest, NaiveCylindricalPathFinderZDirTest)
{
    // set parameters:
    std::map<std::string, real> params;
    params["pfCylNumSteps"] = 10;
    params["pfCylStepLength"] = 0.01;
    params["pfCylRad"] = 1;

    // set direction vector and centre point:
    gmx::RVec dirVec(0.0, 0.0, 2.0*std::sqrt(3.0));
    gmx::RVec centrePoint(0.0, 0.0, 0.0);

    // create path finding module and create path:
    NaiveCylindricalPathFinder pfm(params, centrePoint, dirVec);
    pfm.findPath();

    // find path:
    MolecularPath cylPath = pfm.getMolecularPath();

    // get path radii and support points:
    std::vector<real> radii = cylPath.pathRadii();
    std::vector<gmx::RVec> points = cylPath.pathPoints();

    // assert correct number of points:
    ASSERT_EQ(2*params["pfCylNumSteps"] + 1, radii.size());

    // assert correct radii:
    for(unsigned int i = 0; i < radii.size(); i++)
    {
        ASSERT_NEAR(params["pfCylRad"],
                    radii[i],
                    std::numeric_limits<real>::epsilon());
    }

    // assert correct location of path points:
    for(unsigned int i = 0; i < radii.size(); i++)
    {
        // x-coordinate:
        ASSERT_NEAR(centrePoint[XX],
                    points[i][XX],
                    std::numeric_limits<real>::epsilon());

        // y-coordinate:
        ASSERT_NEAR(centrePoint[YY],
                    points[i][YY],
                    std::numeric_limits<real>::epsilon());

        // z-coordinate:
        ASSERT_NEAR(centrePoint[ZZ] + (i - params["pfCylNumSteps"])*params["pfCylStepLength"],
                    points[i][ZZ],
                    std::numeric_limits<real>::epsilon());
    }
}

