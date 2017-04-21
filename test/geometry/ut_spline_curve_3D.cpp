#include <vector>
#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "geometry/spline_curve_3D.hpp"


/*
 * Test fixture for the one dimensional spline curve.
 */
class SplineCurve3DTest : public ::testing::Test
{
    public:
  
};


/*
 * Uses a simple linear interpolating spline to test whether naive evaluation 
 * works correctly. Evaluation points are chosen to be the control points and
 * the interval midpoints.
 */

TEST_F(SplineCurve3DTest, SplineCurve3DLinearNaiveTest)
{
    // floating point comparison threshold:
    real eps = std::numeric_limits<real>::epsilon();

    // linear spline:
    int degree = 1;

    // evaluate spline itself:
    unsigned int derivOrder = 0;

    // define data points for linear relation:
    std::vector<real> t = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<gmx::RVec> f = {gmx::RVec(-2.0,  2.0,  2.0),
                                gmx::RVec(-1.0,  1.0,  2.5),
                                gmx::RVec( 0.0,  0.0,  3.0),
                                gmx::RVec( 1.0, -1.0,  3.5),
                                gmx::RVec( 2.0, -2.0,  4.0)};

    // create appropriate knot vector for linear interpolation:
    std::vector<real> knots;
    knots.push_back(t.front());
    for(unsigned int i = 0; i < t.size(); i++)
    {
        knots.push_back(t[i]);
    }
    knots.push_back(t.back());
    
    // create corresponding spline curve:
    SplineCurve3D SplC(degree, knots, f);

    // check if spline is evaluates to control points at original data points:
    for(unsigned int i = 0; i < t.size(); i++)
    {
        gmx::RVec value = SplC.evaluate(t[i], derivOrder, eSplineEvalNaive);
        ASSERT_NEAR(f[i][0], value[0], eps);
        ASSERT_NEAR(f[i][1], value[1], eps);
        ASSERT_NEAR(f[i][2], value[2], eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(unsigned int i = 0; i < t.size() - 1; i++)
    {
        real midpoint = (t[i] + t[i+1])/2.0; 
        gmx::RVec value = SplC.evaluate(midpoint, derivOrder, eSplineEvalNaive);
        ASSERT_NEAR((f[i][0] + f[i+1][0])/2.0, value[0], eps);
        ASSERT_NEAR((f[i][1] + f[i+1][1])/2.0, value[1], eps);
        ASSERT_NEAR((f[i][2] + f[i+1][2])/2.0, value[2], eps);
    }
}


/*
 * Uses a simple linear interpolating spline to test whether de Boor evaluation
 * works correctly. Evaluation points are chosen to be the control points and
 * the interval midpoints.
 */

TEST_F(SplineCurve3DTest, SplineCurve3DLinearDeBoorTest)
{
    // floating point comparison threshold:
    real eps = std::numeric_limits<real>::epsilon();

    // linear spline:
    int degree = 1;

    // evaluate spline itself:
    unsigned int derivOrder = 0;

    // define data points for linear relation:
    std::vector<real> t = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<gmx::RVec> f = {gmx::RVec(-2.0,  2.0,  2.0),
                                gmx::RVec(-1.0,  1.0,  2.5),
                                gmx::RVec( 0.0,  0.0,  3.0),
                                gmx::RVec( 1.0, -1.0,  3.5),
                                gmx::RVec( 2.0, -2.0,  4.0)};

    // create appropriate knot vector for linear interpolation:
    std::vector<real> knots;
    knots.push_back(t.front());
    for(unsigned int i = 0; i < t.size(); i++)
    {
        knots.push_back(t[i]);
    }
    knots.push_back(t.back());
    
    // create corresponding spline curve:
    SplineCurve3D SplC(degree, knots, f);

    // check if spline is evaluates to control points at original data points:
    for(unsigned int i = 0; i < t.size(); i++)
    {
        gmx::RVec value = SplC.evaluate(t[i], derivOrder, eSplineEvalNaive);
        ASSERT_NEAR(f[i][0], value[0], eps);
        ASSERT_NEAR(f[i][1], value[1], eps);
        ASSERT_NEAR(f[i][2], value[2], eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(unsigned int i = 0; i < t.size() - 1; i++)
    {
        real midpoint = (t[i] + t[i+1])/2.0; 
        gmx::RVec value = SplC.evaluate(midpoint, derivOrder, eSplineEvalNaive);
        ASSERT_NEAR((f[i][0] + f[i+1][0])/2.0, value[0], eps);
        ASSERT_NEAR((f[i][1] + f[i+1][1])/2.0, value[1], eps);
        ASSERT_NEAR((f[i][2] + f[i+1][2])/2.0, value[2], eps);
    }
}


/*
 * Simple test for whether the first and second derivative of the spline curve
 * are evaluated correctly on a linear spline curve.
 */
TEST_F(SplineCurve3DTest, SplineCurve3DDerivativeTest)
{
    // floating point comparison threshold:
    real eps = std::numeric_limits<real>::epsilon();

    // linear spline:
    int degree = 1;

    // define data points for linear relation:
    std::vector<real> t = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<gmx::RVec> f = {gmx::RVec(-2.0,  2.0,  2.0),
                                gmx::RVec(-1.0,  1.0,  2.5),
                                gmx::RVec( 0.0,  0.0,  3.0),
                                gmx::RVec( 1.0, -1.0,  3.5),
                                gmx::RVec( 2.0, -2.0,  4.0)};

    // create appropriate knot vector for linear interpolation:
    std::vector<real> knots;
    knots.push_back(t.front());
    for(unsigned int i = 0; i < t.size(); i++)
    {
        knots.push_back(t[i]);
    }
    knots.push_back(t.back());
    
    // create corresponding spline curve:
    SplineCurve3D SplC(degree, knots, f);

    // test first derivative:
    unsigned int derivOrder = 1;

    // check if spline is evaluates to control points at original data points:
    for(unsigned int i = 0; i < t.size(); i++)
    {
        gmx::RVec value = SplC.evaluate(t[i], derivOrder, eSplineEvalNaive);
        ASSERT_NEAR( 1.0, value[0], eps);
        ASSERT_NEAR(-1.0, value[1], eps);
        ASSERT_NEAR( 0.5, value[2], eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(unsigned int i = 0; i < t.size() - 1; i++)
    {
        real midpoint = (t[i] + t[i+1])/2.0; 
        gmx::RVec value = SplC.evaluate(midpoint, derivOrder, eSplineEvalNaive);
        ASSERT_NEAR( 1.0, value[0], eps);
        ASSERT_NEAR(-1.0, value[1], eps);
        ASSERT_NEAR( 0.5, value[2], eps);
    }

    // test second derivative:
    derivOrder = 2;

    // check if spline is evaluates to control points at original data points:
    for(unsigned int i = 0; i < t.size(); i++)
    {
        gmx::RVec value = SplC.evaluate(t[i], derivOrder, eSplineEvalNaive);
        ASSERT_NEAR(0.0, value[0], eps);
        ASSERT_NEAR(0.0, value[1], eps);
        ASSERT_NEAR(0.0, value[2], eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(unsigned int i = 0; i < t.size() - 1; i++)
    {
        real midpoint = (t[i] + t[i+1])/2.0; 
        gmx::RVec value = SplC.evaluate(midpoint, derivOrder, eSplineEvalNaive);
        ASSERT_NEAR(0.0, value[0], eps);
        ASSERT_NEAR(0.0, value[1], eps);
        ASSERT_NEAR(0.0, value[2], eps);
    }
}


/*
 * Test for evaluation of spline outside the knot vector range. In this case, 
 * linear extrapolation is used and this is checked for the curve value as well
 * as its first and second derivative.
 */
TEST_F(SplineCurve3DTest, SplineCurve3DExtrapolationTest)
{
    // floating point comparison threshold:
    real eps = std::numeric_limits<real>::epsilon();

    // linear spline:
    int degree = 1;

    // define data points for linear relation:
    std::vector<real> t = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<gmx::RVec> f = {gmx::RVec(-2.0,  2.0,  2.0),
                                gmx::RVec(-1.0,  1.0,  2.5),
                                gmx::RVec( 0.0,  0.0,  3.0),
                                gmx::RVec( 1.0, -1.0,  3.5),
                                gmx::RVec( 2.0, -2.0,  4.0)};

    // create appropriate knot vector for linear interpolation:
    std::vector<real> knots;
    knots.push_back(t.front());
    for(unsigned int i = 0; i < t.size(); i++)
    {
        knots.push_back(t[i]);
    }
    knots.push_back(t.back());
    
    // create corresponding spline curve:
    SplineCurve3D SplC(degree, knots, f);

    // check evaluation below data range:
    real evalPoint = -4.0;
    gmx::RVec value = SplC.evaluate(evalPoint, 0, eSplineEvalDeBoor);
    ASSERT_NEAR(-4.0, value[0], eps);
    ASSERT_NEAR( 4.0, value[1], eps);
    ASSERT_NEAR( 1.0, value[2], eps);

    gmx::RVec frstDeriv = SplC.evaluate(evalPoint, 1, eSplineEvalDeBoor);
    ASSERT_NEAR( 1.0, frstDeriv[0], eps);
    ASSERT_NEAR(-1.0, frstDeriv[1], eps);
    ASSERT_NEAR( 0.5, frstDeriv[2], eps);

    gmx::RVec scndDeriv = SplC.evaluate(evalPoint, 2, eSplineEvalDeBoor);
    ASSERT_NEAR(0.0, scndDeriv[0], eps);
    ASSERT_NEAR(0.0, scndDeriv[1], eps);
    ASSERT_NEAR(0.0, scndDeriv[2], eps);

    // check evaluation above data range:
    evalPoint = 4.0;
    value = SplC.evaluate(evalPoint, 0, eSplineEvalDeBoor);
    ASSERT_NEAR( 4.0, value[0], eps);
    ASSERT_NEAR(-4.0, value[1], eps);
    ASSERT_NEAR( 5.0, value[2], eps);

    frstDeriv = SplC.evaluate(evalPoint, 1, eSplineEvalDeBoor);
    ASSERT_NEAR( 1.0, frstDeriv[0], eps);
    ASSERT_NEAR(-1.0, frstDeriv[1], eps);
    ASSERT_NEAR( 0.5, frstDeriv[2], eps);

    scndDeriv = SplC.evaluate(evalPoint, 2, eSplineEvalDeBoor);
    ASSERT_NEAR(0.0, scndDeriv[0], eps);
    ASSERT_NEAR(0.0, scndDeriv[1], eps);
    ASSERT_NEAR(0.0, scndDeriv[2], eps); 
}


/*
 * Checks curve length calculation on linear curve:
 */
TEST_F(SplineCurve3DTest, SplineCurve3DLengthTest)
{
    // floating point comparison threshold:
    real eps = 1e-3;

    // linear spline:
    int degree = 1;

    // define data points for linear relation:
    std::vector<real> t = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<gmx::RVec> f = {gmx::RVec(-2.0,  2.0,  2.0),
                                gmx::RVec(-1.0,  1.0,  2.5),
                                gmx::RVec( 0.0,  0.0,  3.0),
                                gmx::RVec( 1.0, -1.0,  3.5),
                                gmx::RVec( 2.0, -2.0,  4.0)};

    // create appropriate knot vector for linear interpolation:
    std::vector<real> knots;
    knots.push_back(t.front());
    for(unsigned int i = 0; i < t.size(); i++)
    {
        knots.push_back(t[i]);
    }
    knots.push_back(t.back());
    
    // create corresponding spline curve:
    SplineCurve3D SplC(degree, knots, f);

    // calculate endpoint distance:
    real dX = std::abs(f.front()[0] - f.back()[0]);
    real dY = std::abs(f.front()[1] - f.back()[1]);
    real dZ = std::abs(f.front()[2] - f.back()[2]);
    real trueLength = std::sqrt(dX*dX + dY*dY + dZ*dZ);

    // calculate curve length internally:
    real length = SplC.length(eps);

    // both should be equal for linear curve:
    ASSERT_NEAR(trueLength, length, eps);
}

