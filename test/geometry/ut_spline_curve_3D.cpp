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


#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>

#include <gtest/gtest.h>

#include "geometry/spline_curve_3D.hpp"
#include "geometry/cubic_spline_interp_3D.hpp"


/*!
 * \brief Test fixture for the three dimensional spline curve.
 */
class SplineCurve3DTest : public ::testing::Test
{
    public:
  
};


/*!
 * Tests linear spline in three dimensions.
 */
TEST_F(SplineCurve3DTest, SplineCurve3DLinearTest)
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
        gmx::RVec value = SplC.evaluate(t[i], derivOrder);
        ASSERT_NEAR(f[i][XX], value[XX], eps);
        ASSERT_NEAR(f[i][YY], value[YY], eps);
        ASSERT_NEAR(f[i][ZZ], value[ZZ], eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(unsigned int i = 0; i < t.size() - 1; i++)
    {
        real midpoint = (t[i] + t[i+1])/2.0; 
        gmx::RVec value = SplC.evaluate(midpoint, derivOrder);
        ASSERT_NEAR((f[i][XX] + f[i+1][XX])/2.0, value[XX], eps);
        ASSERT_NEAR((f[i][YY] + f[i+1][YY])/2.0, value[YY], eps);
        ASSERT_NEAR((f[i][ZZ] + f[i+1][ZZ])/2.0, value[ZZ], eps);
    }
}


/*!
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
        gmx::RVec value = SplC.evaluate(t[i], derivOrder);
        ASSERT_NEAR( 1.0, value[XX], eps);
        ASSERT_NEAR(-1.0, value[YY], eps);
        ASSERT_NEAR( 0.5, value[ZZ], eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(unsigned int i = 0; i < t.size() - 1; i++)
    {
        real midpoint = (t[i] + t[i+1])/2.0; 
        gmx::RVec value = SplC.evaluate(midpoint, derivOrder);
        ASSERT_NEAR( 1.0, value[XX], eps);
        ASSERT_NEAR(-1.0, value[YY], eps);
        ASSERT_NEAR( 0.5, value[ZZ], eps);
    }

    // test second derivative:
    derivOrder = 2;

    // check if spline is evaluates to control points at original data points:
    for(unsigned int i = 0; i < t.size(); i++)
    {
        gmx::RVec value = SplC.evaluate(t[i], derivOrder);
        ASSERT_NEAR(0.0, value[XX], eps);
        ASSERT_NEAR(0.0, value[YY], eps);
        ASSERT_NEAR(0.0, value[ZZ], eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(unsigned int i = 0; i < t.size() - 1; i++)
    {
        real midpoint = (t[i] + t[i+1])/2.0; 
        gmx::RVec value = SplC.evaluate(midpoint, derivOrder);
        ASSERT_NEAR(0.0, value[XX], eps);
        ASSERT_NEAR(0.0, value[YY], eps);
        ASSERT_NEAR(0.0, value[ZZ], eps);
    }
}


/*!
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
    gmx::RVec value = SplC.evaluate(evalPoint, 0);
    ASSERT_NEAR(-4.0, value[XX], eps);
    ASSERT_NEAR( 4.0, value[YY], eps);
    ASSERT_NEAR( 1.0, value[ZZ], eps);

    gmx::RVec frstDeriv = SplC.evaluate(evalPoint, 1);
    ASSERT_NEAR( 1.0, frstDeriv[XX], eps);
    ASSERT_NEAR(-1.0, frstDeriv[YY], eps);
    ASSERT_NEAR( 0.5, frstDeriv[ZZ], eps);

    gmx::RVec scndDeriv = SplC.evaluate(evalPoint, 2);
    ASSERT_NEAR(0.0, scndDeriv[XX], eps);
    ASSERT_NEAR(0.0, scndDeriv[YY], eps);
    ASSERT_NEAR(0.0, scndDeriv[ZZ], eps);

    // check evaluation above data range:
    evalPoint = 4.0;
    value = SplC.evaluate(evalPoint, 0);
    ASSERT_NEAR( 4.0, value[XX], eps);
    ASSERT_NEAR(-4.0, value[YY], eps);
    ASSERT_NEAR( 5.0, value[ZZ], eps);

    frstDeriv = SplC.evaluate(evalPoint, 1);
    ASSERT_NEAR( 1.0, frstDeriv[XX], eps);
    ASSERT_NEAR(-1.0, frstDeriv[YY], eps);
    ASSERT_NEAR( 0.5, frstDeriv[ZZ], eps);

    scndDeriv = SplC.evaluate(evalPoint, 2);
    ASSERT_NEAR(0.0, scndDeriv[XX], eps);
    ASSERT_NEAR(0.0, scndDeriv[YY], eps);
    ASSERT_NEAR(0.0, scndDeriv[ZZ], eps); 
}


/*!
 * Tests that curve length is determined correctly by creating an interpolating
 * spline on a point set sampled from a helix and comparing the length to the
 * value of the analytical expression.
 */
TEST_F(SplineCurve3DTest, SplineCurve3DLengthTest)
{
    // floating point comparison threshold:
    real eps = 1e-5;

    // define helix parameters:
    const real PI = std::acos(-1.0);
    real tStart = 0.0;
    real tEnd = 2.0*PI;
    real a = 1.573;
    real b = 0.875/(tEnd - tStart);

    // create a point set describing a helix:
    size_t nParams = 1e2;
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
        gmx::RVec splVal = SplC.evaluate(evalPoint, 0);
        gmx::RVec anaVal(a*std::cos(evalPoint), 
                         a*std::sin(evalPoint), 
                         b*evalPoint);

        // these should be the same:
        ASSERT_NEAR(anaVal[XX], splVal[XX], eps);
        ASSERT_NEAR(anaVal[YY], splVal[YY], eps);
        ASSERT_NEAR(anaVal[ZZ], splVal[ZZ], eps);
    }

    // check if length matches analytical expression:
    ASSERT_NEAR((tEnd - tStart)*std::sqrt(a*a + b*b),
                SplC.length(), 
                eps);
}


/*!
 * Tests differential properties of 3D spline curve.
 */
TEST_F(SplineCurve3DTest, SplineCurve3DDifferentialPropertiesTest)
{
    // floating point comparison threshold:
    real eps = 10*std::numeric_limits<real>::epsilon();
    
    // define helix parameters:
    const real PI = std::acos(-1.0);
    real tStart = 0.0;
    real tEnd = 2.0*PI;
    real a = 5.571;
    real b = 1.223/(tEnd - tStart);

    // create a point set describing a helix:
    size_t nParams = 100;
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
        gmx::RVec splVal = SplC.evaluate(evalPoint, 0);
        gmx::RVec anaVal(a*std::cos(evalPoint), 
                         a*std::sin(evalPoint), 
                         b*evalPoint);

        // these should be the same:
        ASSERT_NEAR(anaVal[XX], splVal[XX], eps);
        ASSERT_NEAR(anaVal[YY], splVal[YY], eps);
        ASSERT_NEAR(anaVal[ZZ], splVal[ZZ], eps);
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


/*!
 * Tests whether the re-parameterisation yields a curve with unit speed. To this
 * end a spline curve is constructed by interpolating a point set sampled from
 * a logarithmic spiral with constant z-velocity, re-parameterisation is 
 * employed and the result velocity is checked to be close to one at several 
 * evaluation points.
 */
TEST_F(SplineCurve3DTest, SplineCurve3DArcLengthReparameterisationTest)
{
    // floating point comparison threshold:
    real eps = 2.0*std::sqrt(std::numeric_limits<real>::epsilon());

    // define helix parameters:
    const real PI = std::acos(-1.0);
    real tStart = -2.0*PI;
    real tEnd = 2.0*PI;
    real a = 2.0;
    real b = 1.0/(tEnd - tStart);
    real c = 0.005;

    // create a point set describing a logarithmically spiral helix:
    size_t nParams = 15;
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


/*!
 * Test for the projection of points in Cartesian coordinates onto a spline 
 * curve. Two cases are considered: a linear spline curve and a spline curve 
 * interpolating points sampled from a planar circle. For the linear spline,
 * an extrapolation case is considered as well. Agreement between numerical 
 * procedure and theoretical expectation is accepted if bot the position along
 * the curve and the distance from the curve agree to within two times the 
 * square root of the machine precision.
 */
TEST_F(SplineCurve3DTest, CartesianToCurvilinearInternalTest)
{

    // floating point comparison threshold:
    // NOTE that sqrt(machine epsilon) is the theoretical best precision that
    // can be achieved!
    real eps = 1.1*std::sqrt(std::numeric_limits<real>::epsilon());

    // linear spline:
    int degree = 1;

    // define data points for linear relation (in xy-plane):
    std::vector<real> t = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<gmx::RVec> f = {gmx::RVec( 0.0,  0.0, -2.0),
                                gmx::RVec( 0.0,  0.0, -1.0),
                                gmx::RVec( 0.0,  0.0,  0.0),
                                gmx::RVec( 0.0,  0.0,  1.0),
                                gmx::RVec( 0.0,  0.0,  2.0)};

    // create appropriate knot vector for linear interpolation:
    std::vector<real> knots;
    knots.push_back(t.front());
    for(unsigned int i = 0; i < t.size(); i++)
    {
        knots.push_back(t[i]);
    }
    knots.push_back(t.back());
    
    // create corresponding spline curve:
    // (linear, so no need to re-parameterise)
    SplineCurve3D Spl(degree, knots, f);

    // check that all original points are found to lie on the curve:
    for(unsigned int i = 0; i < f.size(); i++)
    {
        // evaluate curvilinear coordinates of given point:
        gmx::RVec curvi = Spl.cartesianToCurvilinear(f.at(i));

        // check identity with analytical solution:
        ASSERT_NEAR(t[i], curvi[SS], eps);
        ASSERT_NEAR(0.0, curvi[RR], eps);    
    }

    // prepare a set of test points:
    std::vector<gmx::RVec> pts = {gmx::RVec(0.0,  0.0, -0.5),
                                  gmx::RVec(0.0,  0.0,  1.257689),
                                  gmx::RVec(1.0,  0.0, -1.0),
                                  gmx::RVec(0.0, -1.0,  1.0),
                                  gmx::RVec(0.0,  0.0,  5.0),
                                  gmx::RVec(2.0,  2.0, -6.3)};
    std::vector<real> sTrue = {-0.5,
                                1.257689,
                               -1.0,
                                1.0,
                                5.0,
                               -6.3};
    std::vector<real> dTrue = { 0.0,
                                0.0,
                                1.0,
                                1.0,
                                0.0,
                                8.0};
    
    // loop over all test points:
    for(unsigned int i = 0; i < pts.size(); i++)
    {   
        // evaluate curvilinear coordinates of given point:
        gmx::RVec curvi = Spl.cartesianToCurvilinear(pts.at(i));

        // check identity with analytical solution:
        ASSERT_NEAR(sTrue[i], curvi[SS], eps);
        ASSERT_NEAR(dTrue[i], curvi[RR], eps);   
    }

    // cubic spline:
    degree = 3;

    // define curve parameters:
    const real PI = std::acos(-1.0);
    real tStart = 0.0;
    real tEnd = 1.0*PI;
    real a = 1.0;

    // create a point set describing a circle:
    size_t nParams = 1000;
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
                                   0.0)); 
    }

    // create spline by interpolation:
    // (this will already be arc length parameterised)
    CubicSplineInterp3D Interp;
    Spl = Interp(params, points, eSplineInterpBoundaryHermite);

    // check that all original points are found to lie on the curve:
    for(unsigned int i = 0; i < points.size(); i++)
    {
        // evaluate curvilinear coordinates of given point:
        gmx::RVec curvi = Spl.cartesianToCurvilinear(points.at(i));

        // check identity with analytical solution:
        ASSERT_NEAR(params[i], curvi[SS], eps);
        ASSERT_NEAR(0.0, curvi[RR], eps);    
    }

    // define a new set of test points:
    real par = PI/2.0;
    pts = {gmx::RVec(a*std::cos(par), a*std::sin(par), 0.0),
           gmx::RVec(a*std::cos(par), a*std::sin(par), -2.5),
           gmx::RVec(2.0*a*std::cos(par), 2.0*a*std::sin(par), 0.0),
           gmx::RVec(0.5*a*std::cos(par), 0.5*a*std::sin(par), 0.0),
           gmx::RVec(a*std::cos(2.0*par), a*std::sin(2.0*par), 0.0)
           };

    // corresponding spline coordinates:
    sTrue = {par, 
             par, 
             par, 
             par, 
             2.0f*par
             };
    dTrue = {0.0f, 
             2.5f*2.5f, 
             a*a, 
             (0.5f*a)*(0.5f*a), 
             0.0f
             };

    // check that test points are evaluated correctly:
    for(unsigned int i = 0; i < pts.size(); i++)
    {
        // evaluate curvilinear coordinates of given point:
        gmx::RVec curvi = Spl.cartesianToCurvilinear(pts.at(i));

        // check identity with analytical solution:
        ASSERT_NEAR(sTrue[i], curvi[SS], eps);
        ASSERT_NEAR(dTrue[i], curvi[RR], eps);
    }   
}


/*!
 * Tests the mapping of points in Cartesian coordinates onto a spline curve
 * for cases where the test points are located beyond the upper and lower
 * endpoint of the spline curve.
 */
TEST_F(SplineCurve3DTest, CartesianToCurvilinearExternalTest)
{

    // floating point comparison threshold:
    // NOTE that sqrt(machine epsilon) is the theoretical best precision that
    // can be achieved!
    real eps = 1.1*std::sqrt(std::numeric_limits<real>::epsilon());

    // define curve parameters:
    const real PI = std::acos(-1.0);
    real tStart = 0.0;
    real tEnd = 0.5*PI;
    real a = 1.0;

    // create a point set describing a circle:
    size_t nParams = 500; // large number of points for low interpolation error
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
                                   0.0)); 
    }

    // create spline by interpolation:
    // (this will already be arc length parameterised)
    CubicSplineInterp3D Interp;
    SplineCurve3D spl = Interp(params, points, eSplineInterpBoundaryHermite);

    // calculate segment length:
    real len = spl.length();

    // check that all original points are found to lie on the curve:
    for(unsigned int i = 0; i < points.size(); i++)
    {
        // evaluate curvilinear coordinates of given point:
        gmx::RVec curvi = spl.cartesianToCurvilinear(points.at(i));

        // check identity with analytical solution:
        ASSERT_NEAR(params[i], curvi[SS], eps);
        ASSERT_NEAR(0.0, curvi[RR], eps);    
    }

    // define a new set of test points:
    // NOTE: the spline curve is a quarter circle in the top right quadrant
    std::vector<gmx::RVec> pts = {
            gmx::RVec(a*std::cos(0.0), a*std::sin(0.0), 0.0),   // lower endpoint
            gmx::RVec(a*std::cos(1.0), a*std::sin(1.0), 0.0),   // upper endpoint
            gmx::RVec(a, -1.0, 0.0),                            // beyond lower endpoint
            gmx::RVec(a, -2.0, 0.0),                            // beyond lower endpoint
            gmx::RVec(a, -3.0, 0.0),                            // beyond lower endpoint
            gmx::RVec(a, -1.0, 1.0),                            // beyond lower endpoint with offset
            gmx::RVec(a, -2.0, 0.5),                            // beyond lower endpoint with offset
            gmx::RVec(a, -3.0, 1.5),                            // beyond lower endpoint with offset
            gmx::RVec(-1.0, a, 0.0),                            // beyond upper endpoint
            gmx::RVec(-2.0, a, 0.0),                            // beyond upper endpoint
            gmx::RVec(-3.0, a, 0.0),                            // beyond upper endpoint
            gmx::RVec(-1.0, a, 1.0),                            // beyond upper endpoint with offset
            gmx::RVec(-2.0, a, 0.5),                            // beyond upper endpoint with offset
            gmx::RVec(-3.0, a, 1.5),                            // beyond upper endpoint with offset
            };

    // corresponding spline coordinates:
    std::vector<real> sTrue = {
            0.0, 
            1.0, 
            -1.0, 
            -2.0, 
            -3.0, 
            -1.0, 
            -2.0, 
            -3.0, 
            1.0f + len, 
            2.0f + len, 
            3.0f + len, 
            1.0f + len, 
            2.0f + len, 
            3.0f + len, 
            };
    std::vector<real> dTrue = {
            0.0, 
            0.0, 
            0.0, 
            0.0, 
            0.0, 
            1.0*1.0, 
            0.5*0.5, 
            1.5*1.5, 
            0.0, 
            0.0, 
            0.0, 
            1.0*1.0, 
            0.5*0.5, 
            1.5*1.5, 
            };

    // check that test points are evaluated correctly:
    for(unsigned int i = 0; i < pts.size(); i++)
    {
        // evaluate curvilinear coordinates of given point:
        gmx::RVec curvi = spl.cartesianToCurvilinear(pts.at(i));

        // check identity with analytical solution:
        ASSERT_NEAR(sTrue[i], curvi[SS], eps);
        ASSERT_NEAR(dTrue[i], curvi[RR], eps);
    }   
}

