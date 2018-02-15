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

#include <gtest/gtest.h>

#include "geometry/spline_curve_1D.hpp"


/*!
 * \brief Test fixture for the one dimensional spline curve.
 */
class SplineCurve1DTest : public ::testing::Test
{
    public:
 
        std::vector<real> prepareKnotVector(
                const std::vector<real> &uniqueKnots, 
                unsigned int degree)
        {
            std::vector<real> knots;
            for(unsigned int i = 0; i < degree; i++)
            {
                knots.push_back(uniqueKnots.front());
            }
            for(unsigned int i = 0; i < uniqueKnots.size(); i++)
            {
                knots.push_back(uniqueKnots[i]);
            }
            for(unsigned int i = 0; i < degree; i++)
            {
                knots.push_back(uniqueKnots.back());
            }

            return knots;
        }
 
};


/*!
 * Simple test case for the interval finding method.
 */
TEST_F(SplineCurve1DTest, SplineCurve1DFindIntervalTest)
{
    // use cubic spline here:
    int degree = 3;

    // set up knot vector:
    std::vector<real> knotVector = {-1.0, -0.5, -0.25, -0.2, 0.2, -0.25, 0.5, 1.0};

    // set up appropriate number of control points (value is irrelevant here):
    std::vector<real> ctrlPoints;
    ctrlPoints.resize(knotVector.size() - degree - 1);

    // create spline curve:
    SplineCurve1D SplC(degree, knotVector, ctrlPoints);

    // note special treatment for final point:
    std::vector<real> evalPoints = {-1.0, -0.75, -0.5, 0.0, 0.15, 0.5, 1.0};
    std::vector<int> refIndices = {0, 0, 1, 3, 3, 6, 6};

    // loop over evaluation points:
    for(unsigned int i = 0; i < evalPoints.size(); i++)
    {
        // check that correct interval is found:
        ASSERT_EQ(refIndices[i], SplC.findInterval(evalPoints[i]));
    }
}


/*!
 * Uses a simple linear interpolating spline to test whether spline evaluation 
 * works correctly. Evaluation points are chosen to be the control points and
 * the interval midpoints.
 */
TEST_F(SplineCurve1DTest, SplineCurve1DLinearTest)
{
    // floating point comparison threshold:
    real eps = std::numeric_limits<real>::epsilon();

    // linear spline:
    int degree = 1;

    // evaluate spline itself:
    unsigned int derivOrder = 0;

    // define data points for linear relation:
    std::vector<real> x = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<real> y = {-2.0, -1.0, 0.0, 1.0, 2.0};

    // create appropriate knot vector for linear interpolation:
    std::vector<real> knots;
    knots.push_back(x.front());
    for(unsigned int i = 0; i < x.size(); i++)
    {
        knots.push_back(x[i]);
    }
    knots.push_back(x.back());

    // create corresponding spline curve:
    SplineCurve1D SplC(degree, knots, y);

    // check if spline is evaluates to control points at original data points:
    for(unsigned int i = 0; i < x.size(); i++)
    {
        ASSERT_NEAR(y[i],
                    SplC.evaluate(x[i], derivOrder), 
                    eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(unsigned int i = 0; i < x.size() - 1; i++)
    {
        real midpoint = (x[i] + x[i+1])/2.0; 
        ASSERT_NEAR((y[i] + y[i+1])/2.0, 
                    SplC.evaluate(midpoint, derivOrder), 
                    eps);
    }
}


/*!
 * Simple test for whether the first and second derivative of the spline curve
 * are evaluated correctly on a linear spline curve.
 */
TEST_F(SplineCurve1DTest, SplineCurve1DDerivativeTest)
{
    // floating point comparison threshold:
    real eps = std::numeric_limits<real>::epsilon();

    // linear spline:
    int degree = 1;
    
    // define data points for linear relation:
    std::vector<real> x = {-2.0, -1.0, 0.0,  1.0,  2.0};
    std::vector<real> y = { 4.0,  2.0, 0.0, -2.0, -4.0};

    // create appropriate knot vector for linear interpolation:
    std::vector<real> knots;
    knots.push_back(x.front());
    for(unsigned int i = 0; i < x.size(); i++)
    {
        knots.push_back(x[i]);
    }
    knots.push_back(x.back());

    // create corresponding spline curve:
    SplineCurve1D SplC(degree, knots, y);

    // check if spline is evaluates to control points at original data points:
    for(unsigned int i = 1; i < x.size() - 1; i++)
    {
        real frstDeriv = SplC.evaluate(x[i], 1);
        real scndDeriv = SplC.evaluate(x[i], 2);
        ASSERT_NEAR(-2.0, frstDeriv, eps);
        ASSERT_NEAR(0.0, scndDeriv, eps);
    }

    // check evaluation ad interval midpoints:
    for(unsigned int i = 0; i < x.size() - 1; i++)
    {
        real midpoint = (x[i] + x[i+1])/2.0; 
        real frstDeriv = SplC.evaluate(midpoint, 1);
        real scndDeriv = SplC.evaluate(midpoint, 2);
        ASSERT_NEAR(-2.0, frstDeriv, eps);
        ASSERT_NEAR(0.0, scndDeriv, eps);
    }
}


/*!
 * Test for evaluation of spline outside the knot vector range. In this case, 
 * linear extrapolation is used and this is checked for the curve value as well
 * as its first and second derivative.
 */
TEST_F(SplineCurve1DTest, SplineCurve1DExtrapolationTest)
{
    // floating point comparison threshold:
    real eps = std::numeric_limits<real>::epsilon();

    // linear spline:
    int degree = 1;
    
    // define data points for linear relation:
    std::vector<real> x = {-2.0, -1.0, 0.0,  1.0,  2.0};
    std::vector<real> y = { 4.0,  2.0, 0.0, -2.0, -4.0};

    // create appropriate knot vector for linear interpolation:
    std::vector<real> knots;
    knots.push_back(x.front());
    for(unsigned int i = 0; i < x.size(); i++)
    {
        knots.push_back(x[i]);
    }
    knots.push_back(x.back());

    // create corresponding spline curve:
    SplineCurve1D SplC(degree, knots, y);

    // check evaluation below data range:
    real evalPoint = -4.0;
    ASSERT_NEAR(4.0,
                SplC.evaluate(evalPoint, 0),
                eps);
    ASSERT_NEAR(0.0,
                SplC.evaluate(evalPoint, 1),
                eps);
    ASSERT_NEAR(0.0,
                SplC.evaluate(evalPoint, 2),
                eps);

    // check evaluation above data range:
    evalPoint = 4.0;
    ASSERT_NEAR(-4.0,
                SplC.evaluate(evalPoint, 0),
                eps);
    ASSERT_NEAR(0.0,
                SplC.evaluate(evalPoint, 1),
                eps);
    ASSERT_NEAR(0.0,
                SplC.evaluate(evalPoint, 2),
                eps); 
}

