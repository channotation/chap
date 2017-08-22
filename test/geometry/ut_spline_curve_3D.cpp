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


/*
 *
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
        ASSERT_NEAR(f[i][0], value[0], eps);
        ASSERT_NEAR(f[i][1], value[1], eps);
        ASSERT_NEAR(f[i][2], value[2], eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(unsigned int i = 0; i < t.size() - 1; i++)
    {
        real midpoint = (t[i] + t[i+1])/2.0; 
        gmx::RVec value = SplC.evaluate(midpoint, derivOrder);
        ASSERT_NEAR((f[i][0] + f[i+1][0])/2.0, value[0], eps);
        ASSERT_NEAR((f[i][1] + f[i+1][1])/2.0, value[1], eps);
        ASSERT_NEAR((f[i][2] + f[i+1][2])/2.0, value[2], eps);
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
        ASSERT_NEAR( 1.0, value[0], eps);
        ASSERT_NEAR(-1.0, value[1], eps);
        ASSERT_NEAR( 0.5, value[2], eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(unsigned int i = 0; i < t.size() - 1; i++)
    {
        real midpoint = (t[i] + t[i+1])/2.0; 
        gmx::RVec value = SplC.evaluate(midpoint, derivOrder);
        ASSERT_NEAR( 1.0, value[0], eps);
        ASSERT_NEAR(-1.0, value[1], eps);
        ASSERT_NEAR( 0.5, value[2], eps);
    }

    // test second derivative:
    derivOrder = 2;

    // check if spline is evaluates to control points at original data points:
    for(unsigned int i = 0; i < t.size(); i++)
    {
        gmx::RVec value = SplC.evaluate(t[i], derivOrder);
        ASSERT_NEAR(0.0, value[0], eps);
        ASSERT_NEAR(0.0, value[1], eps);
        ASSERT_NEAR(0.0, value[2], eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(unsigned int i = 0; i < t.size() - 1; i++)
    {
        real midpoint = (t[i] + t[i+1])/2.0; 
        gmx::RVec value = SplC.evaluate(midpoint, derivOrder);
        ASSERT_NEAR(0.0, value[0], eps);
        ASSERT_NEAR(0.0, value[1], eps);
        ASSERT_NEAR(0.0, value[2], eps);
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
    ASSERT_NEAR(-4.0, value[0], eps);
    ASSERT_NEAR( 4.0, value[1], eps);
    ASSERT_NEAR( 1.0, value[2], eps);

    gmx::RVec frstDeriv = SplC.evaluate(evalPoint, 1);
    ASSERT_NEAR( 1.0, frstDeriv[0], eps);
    ASSERT_NEAR(-1.0, frstDeriv[1], eps);
    ASSERT_NEAR( 0.5, frstDeriv[2], eps);

    gmx::RVec scndDeriv = SplC.evaluate(evalPoint, 2);
    ASSERT_NEAR(0.0, scndDeriv[0], eps);
    ASSERT_NEAR(0.0, scndDeriv[1], eps);
    ASSERT_NEAR(0.0, scndDeriv[2], eps);

    // check evaluation above data range:
    evalPoint = 4.0;
    value = SplC.evaluate(evalPoint, 0);
    ASSERT_NEAR( 4.0, value[0], eps);
    ASSERT_NEAR(-4.0, value[1], eps);
    ASSERT_NEAR( 5.0, value[2], eps);

    frstDeriv = SplC.evaluate(evalPoint, 1);
    ASSERT_NEAR( 1.0, frstDeriv[0], eps);
    ASSERT_NEAR(-1.0, frstDeriv[1], eps);
    ASSERT_NEAR( 0.5, frstDeriv[2], eps);

    scndDeriv = SplC.evaluate(evalPoint, 2);
    ASSERT_NEAR(0.0, scndDeriv[0], eps);
    ASSERT_NEAR(0.0, scndDeriv[1], eps);
    ASSERT_NEAR(0.0, scndDeriv[2], eps); 
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


/*!
 * Tests whether the reparameterisation yields a curve with unit speed. To this
 * end a spline curve is constructed by interpolating a point set sampled from
 * a logarithmic spiral with constant z-velocity, reparameterisation is 
 * emploeyed and the result velocity is checked to be close to one at several 
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
 * Test for the projection of points in cartesian coordinatas onto a spline 
 * curve. Two cases are considered: a linear spline curve and a spline curve 
 * interpolating points sampled from a planar circle. For the linear spline,
 * an extrapolation case is considered as well. Agreement between numerical 
 * procedure and theoretical expectation is accepted if bot the position along
 * the curve and the distance from the curve agree to within two times the 
 * square root of the machine precision.
 */
TEST_F(SplineCurve3DTest, CartesianToCurvilinearTest)
{

    // floating point comparison threshold:
    // NOTE that sqrt(machine epsilon) is the theoretical best precision that
    // can be achieved!
    real eps = 2.0*std::sqrt(std::numeric_limits<real>::epsilon());

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
    // (linear, so no need to reparameterise)
    SplineCurve3D Spl(degree, knots, f);

    // check that all original points are found to lie on the curve:
    for(unsigned int i = 0; i < f.size(); i++)
    {
        // evaluate curvilinear coordinates of given point:
        gmx::RVec curvi = Spl.cartesianToCurvilinear(f.at(i),
                                                     -2.1,
                                                     2.1);

        // check identity with analytical solution:
        ASSERT_NEAR(t[i], curvi[0], eps);
        ASSERT_NEAR(0.0, curvi[1], eps);    
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
        gmx::RVec curvi = Spl.cartesianToCurvilinear(pts.at(i),
                                                     -10,
                                                     10);

        // check identity with analytical solution:
        ASSERT_NEAR(sTrue[i], curvi[0], eps);
        ASSERT_NEAR(dTrue[i], curvi[1], eps);   
    }

    // cubic spline:
    degree = 3;

    // define curve parameters:
    const real PI = std::acos(-1.0);
    real tStart = 0.0;
    real tEnd = 1.0*PI;
    real a = 1.0;

    // create a point set describing a circle:
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
        gmx::RVec curvi = Spl.cartesianToCurvilinear(points.at(i),
                                                     tStart,
                                                     tEnd);

        std::cerr<<"i = "<<i<<"  "
                 <<"params[i] = "<<params[i]<<"  "
                 <<"curvi[0] = "<<curvi[0]<<"  "
                 <<"eps = "<<params[i] - curvi[0]<<"  "
                 <<"delta = "<<curvi[1]<<"  "
                 <<std::endl;

        // check identity with analytical solution:
      ASSERT_NEAR(params[i], curvi[0], eps);
      ASSERT_NEAR(0.0, curvi[1], eps);    
    }

    // define a new set of test points:
    real par = PI/2.0;
    pts = {gmx::RVec(a*std::cos(par), a*std::sin(par), 0.0),
           gmx::RVec(a*std::cos(par), a*std::sin(par), -2.5),
           gmx::RVec(2.0*a*std::cos(par), 2.0*a*std::sin(par), 0.0),
           gmx::RVec(0.5*a*std::cos(par), 0.5*a*std::sin(par), 0.0),
           gmx::RVec(a*std::cos(2.0*par), a*std::sin(2.0*par), 0.0)};

    // corresponding spline coordinates:
    sTrue = {par, 
             par, 
             par, 
             par, 
             2.0f*par};
    dTrue = {0.0f, 
             2.5f*2.5f, 
             a*a, 
             (0.5f*a)*(0.5f*a), 
             0.0f};

    // check that test points are evaluated correctly:
    for(unsigned int i = 0; i < pts.size(); i++)
    {
        // evaluate curvilinear coordinates of given point:
        gmx::RVec curvi = Spl.cartesianToCurvilinear(pts.at(i),
                                                     0.5*par,
                                                     2.5*par);

        // check identity with analytical solution:
        // FIXME: this test still fails
        ASSERT_NEAR(sTrue[i], curvi[0], eps);
        ASSERT_NEAR(dTrue[i], curvi[1], eps);
    }   
}

