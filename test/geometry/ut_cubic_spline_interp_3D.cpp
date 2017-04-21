#include <gtest/gtest.h>

#include "geometry/cubic_spline_interp_3D.hpp"


/*
 * Test fixture for cubic spline interpolation in 3D.
 */
class CubicSplineInterp3DTest : public ::testing::Test
{


};


/*
 * Tests interpolation algorithm on a problem which should yield a linear 
 * polynomial. Correct evaluation is checked at the support points and the 
 * interval midpoints. Uses Hermite boundary conditions.
 */
TEST_F(CubicSplineInterp3DTest, CubicSplineInterpHermiteLinearTest)
{
    // set tolerance threshold for numerical comparison:
    real eps = std::numeric_limits<real>::epsilon();

    // create a set of points on a line:
    std::vector<gmx::RVec> points = {gmx::RVec( 3.0,  2.0,  0.5),
                                     gmx::RVec( 3.0,  1.0,  0.5),
                                     gmx::RVec( 3.0,  0.0,  0.5),
                                     gmx::RVec( 3.0, -1.0,  0.5),
                                     gmx::RVec( 3.0, -2.0,  0.5)};

    // instantiate interpolation functor:
    CubicSplineInterp3D Interp;

    // find interpolating spline curve:
    SplineCurve3D Spl = Interp(points, eSplineInterpBoundaryHermite);

    // check that spline curve goes through support points:
    for(unsigned int i = 0; i < points.size(); i++)
    {
        real evalPoint = 1.0*i;
        gmx::RVec val = Spl.evaluate(evalPoint, 0, eSplineEvalDeBoor);

        ASSERT_NEAR(points[i][0], val[0], eps);
        ASSERT_NEAR(points[i][1], val[1], eps);
        ASSERT_NEAR(points[i][2], val[2], eps);
    }

    // check interpolation at interval midpoints:
    for(unsigned int i = 0; i < points.size() - 1; i++)
    {
        real evalPoint = 1.0*i + 0.5;
        gmx::RVec val = Spl.evaluate(evalPoint, 0, eSplineEvalDeBoor);

        ASSERT_NEAR((points[i][0] + points[i+1][0])/2.0, val[0], eps);
        ASSERT_NEAR((points[i][1] + points[i+1][1])/2.0, val[1], eps);
        ASSERT_NEAR((points[i][2] + points[i+1][2])/2.0, val[2], eps);
    }
}

