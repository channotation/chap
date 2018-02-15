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

#include "geometry/cubic_spline_interp_1D.hpp"


/*!
 * \brief Test fixture for cubic spline interpolation in 1D.
 */
class CubicSplineInterp1DTest : public ::testing::Test
{


};


/*!
 * Tests the interpolation algorithm on a problem which should yield a linear
 * polynomial. Correct evaluation is checked at the support points and the 
 * interval midpoints. Uses Hermite boundary conditions.
 */
TEST_F(CubicSplineInterp1DTest, CubicSplineInterpHermiteLinearTest)
{
    // will evaluate spline rather than derivatives:
    unsigned int derivOrder = 0;

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
        real val = Spl.evaluate(x.at(i), derivOrder);
        ASSERT_NEAR(f.at(i), val, std::numeric_limits<real>::epsilon());
    }

    // check that spline interpolates correctly at interval midpoints:
    for(unsigned int i = 0; i < x.size() - 1; i++)
    {
        real evalPoint = (x.at(i) + x.at(i+1))/2.0;  
        real val = Spl.evaluate(evalPoint, derivOrder);
        ASSERT_NEAR((f.at(i) + f.at(i+1))/2.0, val, std::numeric_limits<real>::epsilon());
    }
}

