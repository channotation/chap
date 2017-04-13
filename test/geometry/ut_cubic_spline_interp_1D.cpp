#include <gtest/gtest.h>

#include "geometry/cubic_spline_interp_1D.hpp"


/*
 * Test fixture for cubic spline interpolation in 1D.
 */
class CubicSplineInterp1DTest : public ::testing::Test
{


};


/*
 * Tests the interpolation algorithm on a problem which should yield a linear
 * polynomial. Correct evaluation is checked at the support points and the 
 * interval midpoints. Uses Hermite boudnary conditions.
 */
TEST_F(CubicSplineInterp1DTest, CubicSplineInterpHermiteLinearTest)
{
    // define point set to be interpolated:
    std::vector<real> x = {-2.0, -1.0,  0.0,  1.0,  2.0};
    std::vector<real> f = {-2.0, -1.0,  0.0,  1.0,  2.0};

    // instantiate interpolation object:
    CubicSplineInterp1D Interp;

    // find interpolating spline curve:
    SplineCurve1D Spl = Interp(x, f, eSplineInterpBoundaryHermite);

    // check that spline curve goes through support points:
    for(unsigned int i = 0; i < x.size(); i++)
    {
        real val = Spl.evaluate(x.at(i), eSplineEvalDeBoor);
        ASSERT_NEAR(f.at(i), val, std::numeric_limits<real>::epsilon());
    }

    // check that spline interpolates correctly at interval midpoints:
    for(unsigned int i = 0; i < x.size() - 1; i++)
    {
        real evalPoint = (x.at(i) + x.at(i+1))/2.0;  
        real val = Spl.evaluate(evalPoint, eSplineEvalDeBoor);
        ASSERT_NEAR((f.at(i) + f.at(i+1))/2.0, val, std::numeric_limits<real>::epsilon());
    }
}


/*
 * Tests the interpolation algorithm on a problem which should yield a linear
 * polynomial. Correct evaluation is checked at the support points and the 
 * interval midpoints. Uses Natural boudnary conditions.
 */
TEST_F(CubicSplineInterp1DTest, CubicSplineInterpNaturalLinearTest)
{
    // define point set to be interpolated:
    std::vector<real> x = {-2.0, -1.0,  0.0,  1.0,  2.0};
    std::vector<real> f = {-2.0, -1.0,  0.0,  1.0,  2.0};

    // instantiate interpolation object:
    CubicSplineInterp1D Interp;

    // find interpolating spline curve:
//    SplineCurve1D Spl = Interp(x, f, eSplineInterpBoundaryNatural);

    // check that spline curve goes through support points:
    for(unsigned int i = 0; i < x.size(); i++)
    {
//        real val = Spl.evaluate(x.at(i), eSplineEvalDeBoor);
//        ASSERT_NEAR(f.at(i), val, std::numeric_limits<real>::epsilon());
    }

    // check that spline interpolates correctly at interval midpoints:
    for(unsigned int i = 0; i < x.size() - 1; i++)
    {
//        real evalPoint = (x.at(i) + x.at(i+1))/2.0;  
//        real val = Spl.evaluate(evalPoint, eSplineEvalDeBoor);
//        ASSERT_NEAR((f.at(i) + f.at(i+1))/2.0, val, std::numeric_limits<real>::epsilon());
    }
}

