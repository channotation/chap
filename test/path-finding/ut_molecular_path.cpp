#include <gtest/gtest.h>

#include <gromacs/math/vec.h>

#include "path-finding/molecular_path.hpp"


/*
 *
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
         */
        MolecularPath makeHourglassPath(gmx::RVec dir,
                                        gmx::RVec centre,
                                        real length,
                                        real minRadius,
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
                real radius = minRadius + 
                    stepLength*stepLength*(numPoints/2 - i)*(numPoints/2 - i);
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


/*
 *
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
    ASSERT_NEAR(length, mpCylindrical.length(), 10*eps);
    

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
    ASSERT_NEAR(length, mpHourglass.length(), 10*eps);


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


    // create a spring-shaped  path:
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


/*
 *
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
    ASSERT_NEAR(radius, mpCylindrical.minRadius().second, 10*eps);


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

    // assert correct minimal radius and location of minimum:
    ASSERT_NEAR(length/2.0, mpHourglass.minRadius().first, eps);
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


/*
 *
 */
TEST_F(MolecularPathTest, MolecularPathVolumeTest)
{
/*
    // create a toroidal path    
    real pathRadius = 1.0;
    real torusRadius = 10.0;
    real zOffset = 0.0;
    int numPoints = 100;
    MolecularPath mp = makeToroidalPath(pathRadius,
                                        torusRadius,
                                        zOffset,
                                        numPoints);
*/
}

