#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>

#include <gtest/gtest.h>

#include "geometry/spline_curve_3D.hpp"
#include "geometry/cubic_spline_interp_3D.hpp"


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
 * Tests that curve length is determined correctly by creating an interpolating
 * spline on a point set sampled from a helix and comparing the length to the
 * value of the analytical expression.
 */
TEST_F(SplineCurve3DTest, SplineCurve3DLengthTest)
{
    // floating point comparison threshold:
    real eps = 1e-5;

    // cubic spline:
    int degree = 3;

    // define helix parameters:
    const real PI = std::acos(-1.0);
    real tStart = 0.0;
    real tEnd = 2.0*PI;
    real a = 1.573;
    real b = 0.875/(tEnd - tStart);

    // create a point set describing a helix:
    int nParams = 1e2;
    real paramStep = (tEnd - tStart) / (nParams - 1);
    std::vector<real> params;
    std::vector<gmx::RVec> points;
    for(unsigned int i = 0; i < nParams; i++)
    {
        // setup parameter vector:
        params.push_back(i*paramStep + tStart);

        // set up curve:
        points.push_back(gmx::RVec(a*std::cos(params.back()),
                                   a*std::sin(params.back()),
                                   b*params.back())); 
    }

    // create spline by interpolation:
    CubicSplineInterp3D Interp;
    SplineCurve3D SplC = Interp(params, points, eSplineInterpBoundaryHermite);

    // check if spline curve approximates true helix well:
    unsigned int nEval = 100;
    for(unsigned int i = 0; i < nEval; i++)
    {
        // calculate evaluation point:
        real evalPoint = i*(tEnd - tStart)/(nEval - 1);

        // evaluate spline and analytical expression at this point:
        gmx::RVec splVal = SplC(evalPoint, 0, eSplineEvalDeBoor);
        gmx::RVec anaVal(a*std::cos(evalPoint), 
                         a*std::sin(evalPoint), 
                         b*evalPoint);

        // these should be the same:
        ASSERT_NEAR(anaVal[0], splVal[0], eps);
        ASSERT_NEAR(anaVal[1], splVal[1], eps);
        ASSERT_NEAR(anaVal[2], splVal[2], eps);
    }

    // check if length mataches analytical expression:
    ASSERT_NEAR((tEnd - tStart)*std::sqrt(a*a + b*b),
                SplC.length(), 
                eps);
}


/*
 * 
 */
TEST_F(SplineCurve3DTest, SplineCurve3DDifferentialPropertiesTest)
{
    // floating point comparison threshold:
    real eps = 10*std::numeric_limits<real>::epsilon();
    
    // cubic spline:
    int degree = 3;

    // define helix parameters:
    const real PI = std::acos(-1.0);
    real tStart = 0.0;
    real tEnd = 2.0*PI;
    real a = 5.571;
    real b = 1.223/(tEnd - tStart);

    // create a point set describing a helix:
    int nParams = 100;
    real paramStep = (tEnd - tStart) / (nParams - 1);
    std::vector<real> params;
    std::vector<gmx::RVec> points;
    for(unsigned int i = 0; i < nParams; i++)
    {
        // setup parameter vector:
        params.push_back(i*paramStep + tStart);

        // set up curve:
        points.push_back(gmx::RVec(a*std::cos(params.back()),
                                   a*std::sin(params.back()),
                                   b*params.back())); 
    }

    // create spline by interpolation:
    CubicSplineInterp3D Interp;
    SplineCurve3D SplC = Interp(params, points, eSplineInterpBoundaryHermite);

    // check if spline curve approximates true helix well:
    unsigned int nEval = 100;
    for(unsigned int i = 0; i < nEval; i++)
    {
        // calculate evaluation point:
        real evalPoint = i*(tEnd - tStart)/(nEval - 1);

        // evaluate spline and analytical expression at this point:
        gmx::RVec splVal = SplC(evalPoint, 0, eSplineEvalDeBoor);
        gmx::RVec anaVal(a*std::cos(evalPoint), 
                         a*std::sin(evalPoint), 
                         b*evalPoint);

        // these should be the same:
        ASSERT_NEAR(anaVal[0], splVal[0], eps);
        ASSERT_NEAR(anaVal[1], splVal[1], eps);
        ASSERT_NEAR(anaVal[2], splVal[2], eps);
    }

    // check if spline curve speed evaluation is correct:
    nEval = 10;
    for(unsigned int i = 0; i < nEval; i++)
    {
        // calculate evaluation point:
        real evalPoint = i*(tEnd - tStart)/(nEval - 1);

        // evaluate spline and analytical expression at this point:
        real splVal = SplC.speed(evalPoint);
        real anaVal = std::sqrt(a*a + b*b);

        // these should be the same:
        ASSERT_NEAR(anaVal, splVal, 1e-2);
    } 
}


/*
 * Tests whether the reparameterisation yields a curve with unit speed. To this
 * end a spline curve is constructed by interpolating a point set sampled from
 * a logarithmic spiral with constant z-velocity, reparameterisation is 
 * emploeyed and the result velocity is checked to be close to one at several 
 * evaluation points.
 */
TEST_F(SplineCurve3DTest, SplineCurve3DArcLengthReparameterisationTest)
{
    // floating point comparison threshold:
    real eps = 1e-3;

    // cubic spline:
    int degree = 3;

    // define helix parameters:
    const real PI = std::acos(-1.0);
    real tStart = -2.0*PI;
    real tEnd = 2.0*PI;
    real a = 2.0;
    real b = 1.0/(tEnd - tStart);
    real c = 0.005;

    // create a point set describing a logarithmically spiral helix:
    int nParams = 15;
    real paramStep = (tEnd - tStart) / (nParams - 1);
    std::vector<real> params;
    std::vector<gmx::RVec> points;
    for(unsigned int i = 0; i < nParams; i++)
    {
        // setup parameter vector:
        params.push_back(i*paramStep + tStart);

        // set up curve:
        points.push_back(gmx::RVec(a*std::exp(c*params.back())*std::cos(params.back()),
                                   a*std::exp(c*params.back())*std::sin(params.back()),
                                   b*params.back())); 
    }

    // create spline by interpolation:
    CubicSplineInterp3D Interp;
    SplineCurve3D SplC = Interp(params, points, eSplineInterpBoundaryHermite);

    // switch to an arc length parameterised curve:
    SplC.arcLengthParam();

    // check that spline curve now has unit speed everywhere::
    int nEval = 100;
    for(int i = 0; i < nEval; i++)
    {
        // calculate evaluation point:
        real evalPoint = i*(tEnd - tStart)/(nEval - 1);

        // evaluate spline and analytical expression at this point:
        real splVal = SplC.speed(evalPoint);
        real anaVal = 1.0;

        // these should be the same:
        ASSERT_NEAR(anaVal, splVal, eps);
    }     
}

