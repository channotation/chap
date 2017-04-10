#include <vector>
#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "geometry/spline_curve_1D.hpp"


/*
 * Test fucture for the one dimensional spline curve.
 */
class SplineCurve1DTest : public ::testing::Test
{
    public:
  
};



/*
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
    for(int i = 0; i < evalPoints.size(); i++)
    {
        // check that correct interval is found:
        ASSERT_EQ(refIndices[i], SplC.findInterval(evalPoints[i]));
    }
}


/*
 * Uses a simple linear interpolating spline to test whether naive evaluation 
 * works correctly. Evaluation points are chosen to be the control points and
 * the interval midpoints.
 */
TEST_F(SplineCurve1DTest, SplineCurve1DLinearNaiveTest)
{
    // floating point comparison threshold:
    real eps = std::numeric_limits<real>::epsilon();

    // linear spline:
    int degree = 1;

    // define data points for linear relation:
    std::vector<real> x = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<real> y = {-2.0, -1.0, 0.0, 1.0, 2.0};

    // create appropriate knot vector for linear interpolation:
    std::vector<real> knots;
    knots.push_back(x.front());
    for(int i = 0; i < x.size(); i++)
    {
        knots.push_back(x[i]);
    }
    knots.push_back(x.back());

    // create corresponding spline curve:
    SplineCurve1D SplC(degree, knots, y);

    // check if spline is evaluates to control points at original data points:
    for(int i = 0; i < x.size(); i++)
    {
        ASSERT_NEAR(y[i],
                    SplC.evaluate(x[i], eSplineEvalNaive), 
                    eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    for(int i = 0; i < x.size() - 1; i++)
    {
        real midpoint = (x[i] + x[i+1])/2.0; 
        ASSERT_NEAR((y[i] + y[i+1])/2.0, 
                    SplC.evaluate(midpoint, eSplineEvalNaive), 
                    eps);
    }
}


/*
 * Uses a simple linear interpolating spline to test whether de Boor evaluation
 * works correctly. Evaluation points are chosen to be the control points and
 * the interval midpoints.
 */
TEST_F(SplineCurve1DTest, SplineCurve1DLinearDeBoorTest)
{
    // floating point comparison threshold:
    real eps = std::numeric_limits<real>::epsilon();

    // linear spline:
    int degree = 1;

    // define data points for linear relation:
    std::vector<real> x = {-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<real> y = {-2.0, -1.0, 0.0, 1.0, 2.0};

    // create appropriate knot vector for linear interpolation:
    std::vector<real> knots;
    knots.push_back(x.front());
    for(int i = 0; i < x.size(); i++)
    {
        knots.push_back(x[i]);
    }
    knots.push_back(x.back());

    // create corresponding spline curve:
    SplineCurve1D SplC(degree, knots, y);

    // check if spline is evaluates to control points at original data points:
    for(int i = 0; i < x.size(); i++)
    {
        ASSERT_NEAR(y[i],
                    SplC.evaluate(x[i], eSplineEvalDeBoor), 
                    eps);
    }

    // check if spline interpolates linearly at interval midpoints:
    
    for(int i = 0; i < x.size() - 1; i++)
    {
        real midpoint = (x[i] + x[i+1])/2.0; 
        ASSERT_NEAR((y[i] + y[i+1])/2.0, 
                    SplC.evaluate(midpoint, eSplineEvalDeBoor), 
                    eps);
    }
}

