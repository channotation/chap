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

#include <gromacs/math/vec.h>

#include "path-finding/molecular_path.hpp"


/*!
 * \brief Test fixture for MolecularPath objects.
 *
 * Provides some functions to set up simple mathematical shapes as paths,
 * including cylinder, hourglass, torus, and spring.
 */
class MolecularPathTest : public ::testing::Test
{

    public:

        // mathematical constants:
        const real PI_ = std::acos(-1.0);


        /*!
         * Auxiliary function that creates a cylindrical path of given length
         * and radius point in a given direction and centred around a given 
         * position. 
         */
        MolecularPath makeCylindricalPath(gmx::RVec dir,
                                          gmx::RVec centre,
                                          real length,
                                          real radius,
                                          int numPoints)
        {
            // ensure direction is normalised:
            unitv(dir, dir);

            // calculate starting point from centre:
            gmx::RVec startingPoint;
            svmul(-length/2.0, dir, startingPoint);
            rvec_add(startingPoint, centre, startingPoint);
    
            // calculate stepping vector:
            gmx::RVec step;
            svmul(length/(numPoints - 1), dir, step);

            // create vectors of support points and radii:
            std::vector<gmx::RVec> pathPoints;
            std::vector<real> pathRadii;
            pathPoints.push_back(startingPoint);
            pathRadii.push_back(radius);
            for(int i = 1; i < numPoints; i++)
            {
                // move to next point on cylinder line:
                gmx::RVec point(pathPoints[i - 1]);
                rvec_add(point, step, point);
                pathPoints.push_back(point);
                
                // also add to vector of radii:
                pathRadii.push_back(radius);
            }

            // instantiate the path object:
            return MolecularPath(pathPoints, pathRadii);
        }


        /*!
         * Auxiliary function that creates an hourglass shaped path pointing in
         * a given direction and centred around a given position. The centre 
         * line of this path is straight (i.e. zero curvature) and the radius
         * becomes minimal at \f$ s = L/2 \f$. This minimum is only unique if
         * the path is created from an odd number of points.
         *
         * Note that the radius of this pore will be proportional to the third
         * power of the pore length. Such a steep increase is unlikely to 
         * happen in ion channel pores, where the radius will be roughly on the
         * same order of magnitude.
         */
        MolecularPath makeHourglassPath(gmx::RVec dir,
                                        gmx::RVec centre,
                                        real length,
                                        real minRadius,
                                        int numPoints)
        {
            // sanity check:
            if( numPoints % 2 == 0 )
            {
                throw std::runtime_error("ERROR: This function only works for odd number of points!");
                std::abort();
            }
    

            // ensure direction is normalised:
            unitv(dir, dir);

            // calculate starting point from centre:
            gmx::RVec startingPoint;
            svmul(-length/2.0, dir, startingPoint);
            rvec_add(startingPoint, centre, startingPoint);
    
            // calculate stepping vector:
            gmx::RVec step;
            svmul(length/(numPoints - 1), dir, step);
            real stepLength = norm(step);

            // create vectors of support points and radii:
            std::vector<gmx::RVec> pathPoints;
            std::vector<real> pathRadii;
            pathPoints.push_back(startingPoint);
            for(int i = 1; i < numPoints; i++)
            {
                // move to next point on cylinder line:
                gmx::RVec point(pathPoints[i - 1]);
                rvec_add(point, step, point);
                pathPoints.push_back(point);
            }


            // build up radii:
            for(int i = 0; i < numPoints; i++)
            {
                real s = stepLength*i - length/2.0;
                real radius = minRadius + 
                    std::abs(std::pow(s, 3));
                pathRadii.push_back(radius);
            }

            // instantiate the path object:
            return MolecularPath(pathPoints, pathRadii);   
        }


        /*!
         * Creates a path of the shaped of half a torus ("u-shape", 
         * "half donut") and of constant radius along this path. The torus lies
         * in the \f$ xy \f$ -plane and can be shifted by a fixed amount along
         * the z-axis. Note that a full torus shape can not be realised, as in
         * this case the interpolaion problem becomes degenerate.
         */
        MolecularPath makeToroidalPath(real pathRadius,
                                       real torusRadius,
                                       real zOffset,
                                       int numPoints)
        {
            // support point step length:
            real ds = PI_/(numPoints - 1);

            // create vectors of support points and radii:
            std::vector<gmx::RVec> pathPoints;
            std::vector<real> pathRadii;
            for(int i = 0; i < numPoints; i++)
            {
                // sample support point on circle in xy-plane:
                gmx::RVec point(torusRadius*std::cos(ds*i),
                                torusRadius*std::sin(ds*i),
                                zOffset);
                pathPoints.push_back(point);
                
                // also add to vector of radii:
                pathRadii.push_back(pathRadius);
            }

            // instantiate the path object:
            return MolecularPath(pathPoints, pathRadii);
        };


        /*!
         * Creates a molecular path with the shape of a spring, i.e. a tube 
         * with circular cross section wrapped around a helix. The radius of 
         * the rube is constant and the centre of the shape can be shifted by
         * a certain offset.
         */
        MolecularPath makeSpringPath(real pathRadius,
                                     real springA,
                                     real springB,
                                     real springParLen,
                                     gmx::RVec offset,
                                     int numPoints)
        {
            // support point step length:
            real ds = springParLen/(numPoints - 1);

            // create vectors of support points and radii:
            std::vector<gmx::RVec> pathPoints;
            std::vector<real> pathRadii;
            for(int i = 0; i < numPoints; i++)
            {
                // sample support point on circle in xy-plane:
                gmx::RVec point(springA*std::cos(ds*i) + offset[XX],
                                springA*std::sin(ds*i) + offset[YY],
                                springB*ds*i + offset[ZZ]);
                pathPoints.push_back(point);
                
                // also add to vector of radii:
                pathRadii.push_back(pathRadius);
            }

            // instantiate the path object:
            return MolecularPath(pathPoints, pathRadii);
        }
};


/*!
 * Test that MolecularPath objects return the correct (arc) length. This is 
 * tested on a cylinder, an hourglass, a half-torus, and a spring and the 
 * tolarance threshold for floating point comparison is taken to be
 * \f$ \sqrt{\epsilon} \f$.
 */
TEST_F(MolecularPathTest, MolecularPathLengthTest)
{
    // get machine epsilon:
    real eps = std::numeric_limits<real>::epsilon();

    // create a cylindrical path:
    gmx::RVec dir(1.0, 5.0, -2.2);
    gmx::RVec centre(-0.4, 1.5, 0.3);
    real length = 4.5;
    real radius = 0.75;
    int numPoints = 10;
    MolecularPath mpCylindrical = makeCylindricalPath(
            dir, 
            centre, 
            length, 
            radius,
            numPoints);

    // assert correct length:
    ASSERT_NEAR(length, mpCylindrical.length(), std::sqrt(eps));
    

    // create an hourglass-shaped path:
    dir = gmx::RVec(0.0, 0.0, 1.0);
    centre = gmx::RVec(0.4, -2.5, -0.1);
    length = 4.5;
    radius = 0.15;
    numPoints = 11;
    MolecularPath mpHourglass = makeHourglassPath(
            dir, 
            centre, 
            length, 
            radius,
            numPoints);

    // assert correct length:
    ASSERT_NEAR(length, mpHourglass.length(), std::sqrt(eps));


    // create a toroidal path:
    real pathRadius = 0.5;
    real torusRadius = 10.0;
    real zOffset = -5.3;
    numPoints = 25;
    MolecularPath mpToroidal = makeToroidalPath(pathRadius,
                                                torusRadius,
                                                zOffset,
                                                numPoints);

    // assert correct length:
    ASSERT_NEAR(PI_*torusRadius, mpToroidal.length(), std::sqrt(eps));


    // create a spring-shaped path:
    pathRadius = 0.015;
    real springA = 1.0;
    real springB = 3.0;
    real springParLen = 4.0*PI_;
    gmx::RVec offset(-2.0, 0.3, 1.5);
    numPoints = 30;
    MolecularPath mpSpring = makeSpringPath(pathRadius,
                                            springA,
                                            springB,
                                            springParLen,
                                            offset,
                                            numPoints);

    // assert correct length:
    real springLength = std::sqrt(springA*springA + springB*springB)*springParLen;
    ASSERT_NEAR(springLength, mpSpring.length(), std::sqrt(eps));
}


/*!
 * This tests checks that the MolecularPath object returns the correct minimum
 * (and arg min) radius. The minimum is is checked for a a cylinder, an 
 * hourglass, a half-torus, and a spring, where the tolerance threshold is 
 * taken to be \f$ \epsilon \f$. The arg min is only tested on the hourglass 
 * path (it is not unique in all other cases) with a tolerance of 
 * \f$ \sqrt{\epsilon}  \f$ . 
 */
TEST_F(MolecularPathTest, MolecularPathMinRadiusTest)
{
    // get machine epsilon:
    real eps = std::numeric_limits<real>::epsilon();

    // create a cylindrical path:
    gmx::RVec dir(1.0, 5.0, -2.2);
    gmx::RVec centre(0.0, 0.0, 0.0);
    real length = 4.5;
    real radius = 0.75;
    int numPoints = 10;
    MolecularPath mpCylindrical = makeCylindricalPath(
            dir, 
            centre, 
            length, 
            radius,
            numPoints);

    // assert correct minimum radius (arg min is undefined for constant radius):
    ASSERT_NEAR(radius, mpCylindrical.minRadius().second, std::sqrt(eps));


    // create an hourglass-shaped path:
    dir = gmx::RVec(-0.6, 0.5, 1.0);
    centre = gmx::RVec(0.4, -2.5, -0.1);
    length = 2.1;
    radius = 0.5;
    numPoints = 25;
    MolecularPath mpHourglass = makeHourglassPath(
            dir, 
            centre, 
            length, 
            radius,
            numPoints);

    // assert correct minimal radius and location of minimum:
    ASSERT_NEAR(length/2.0, mpHourglass.minRadius().first, 2*std::sqrt(eps));
    ASSERT_NEAR(radius, mpHourglass.minRadius().second, eps);

    // create a toroidal path    
    real pathRadius = 0.5;
    real torusRadius = 10.0;
    real zOffset = -5.3;
    numPoints = 25;
    MolecularPath mpToroidal = makeToroidalPath(pathRadius,
                                                torusRadius,
                                                zOffset,
                                                numPoints);

    // assert correct minimum radius (argmin is undefined for constant radius):
    ASSERT_NEAR(pathRadius, mpToroidal.minRadius().second, std::sqrt(eps));


    // create a spring-shaped  path:
    pathRadius = 0.015;
    real springA = 1.0;
    real springB = 3.0;
    real springParLen = 4.0*PI_;
    gmx::RVec offset(-2.0, 0.3, 1.5);
    numPoints = 25;
    MolecularPath mpSpring = makeSpringPath(pathRadius,
                                            springA,
                                            springB,
                                            springParLen,
                                            offset,
                                            numPoints);

    // assert correct minimum radius (arg min is undefined in this case):
    ASSERT_NEAR(pathRadius, mpSpring.minRadius().second, eps);
}


/*!
 * This test checks that the MolecularPath object calculates the correct volume
 * for a cylinder, hourglass, half-torus, and spring. The tolerance threshold 
 * is \f$ \sqrt{\epsilon} \f$ in all cases.
 */
TEST_F(MolecularPathTest, MolecularPathVolumeTest)
{
    // get machine epsilon:
    real eps = std::numeric_limits<real>::epsilon();

    // create a cylindrical path:
    gmx::RVec dir(1.0, -2.0, 0.5);
    gmx::RVec centre(15.3, -25.0, 10.0);
    real length = 50.0;
    real radius = 0.05;
    int numPoints = 10;
    MolecularPath mpCylindrical = makeCylindricalPath(
            dir, 
            centre, 
            length, 
            radius,
            numPoints);

    // assert correct volume of cylinder:
    ASSERT_NEAR(PI_*radius*radius*length, 
                mpCylindrical.volume(), 
                std::sqrt(eps));


    // create an hourglass-shaped path:
    // (note: if length too large, radius varies over many orders of magnitude
    // so that chord-length approximation in re-parameterisation is poor, 
    // leading to an error in volume estimation, but this should not be 
    // relevant in channels, where chord length should be a good approximation)
    dir = gmx::RVec(0.2, -3.0, 1.0);
    centre = gmx::RVec(-3.3, 4.0, std::sqrt(2.0));
    length = 2.0; 
    radius = 0.015;
    numPoints = 25;
    MolecularPath mpHourglass = makeHourglassPath(
            dir, 
            centre, 
            length, 
            radius,
            numPoints);

    // assert correct volume:
    real volumeHourglass = std::pow(length, 7)/896.0 + 
                           std::pow(length, 4)/32.0*radius +
                           radius*radius*length/2.0;
    volumeHourglass *= PI_*2.0;
    ASSERT_NEAR(volumeHourglass,
                mpHourglass.volume(), 
                std::sqrt(eps));


    // create a toroidal path    
    real pathRadius = 0.5;
    real torusRadius = 10.0;
    real zOffset = -5.3;
    numPoints = 25;
    MolecularPath mpToroidal = makeToroidalPath(pathRadius,
                                                torusRadius,
                                                zOffset,
                                                numPoints);

    // assert correct volume of torus:
    ASSERT_NEAR(PI_*PI_*pathRadius*pathRadius*torusRadius, 
                mpToroidal.volume(), 
                std::sqrt(eps));

    // create a spring-shaped path:
    pathRadius = 0.015;
    real springA = 1.0;
    real springB = 3.0;
    real springParLen = 4.0*PI_;
    gmx::RVec offset(-2.0, 0.3, 1.5);
    numPoints = 25;
    MolecularPath mpSpring = makeSpringPath(pathRadius,
                                            springA,
                                            springB,
                                            springParLen,
                                            offset,
                                            numPoints);

    // assert correct volume:
    ASSERT_NEAR(PI_*PI_*pathRadius*pathRadius*springParLen, 
                mpSpring.volume(), 
                std::sqrt(eps));                
}

