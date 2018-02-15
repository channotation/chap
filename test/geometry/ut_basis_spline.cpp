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
#include <iostream>
#include <limits>
#include <algorithm>

#include <gtest/gtest.h>

#include "geometry/basis_spline.hpp"


/*!
 * \brief Test fixture for testing the BasisSpline class.
 *
 * Reference values have been computed with the bs() method of the spline 
 * library in R. Values of knots reflect the example discussed in "A very short 
 * note on B-splines" by Samiran Sinha. Evalution points are selected so that 
 * both values at and between the original knots are probed.
 */
class BasisSplineTest : public ::testing::Test
{
				protected:

        // knot vector:
        std::vector<real> knotVector_ = {-4, -0.5, 0.0, 0.5, 4};
         
        // evaluation points:
        std::vector<real> evalPoints_ = {-4.0, -2.5, 0.0, 0.5, -1.0, std::sqrt(2.0), 4.0};
};


/*!
 * Tests that the basis splines over a knot vector form a partition of unity,
 * i.e. that sum_i=1^n B_i,k(x) = 1. This is done for basis splines up to
 * degree 5. Test is considered passed if the above sum is within 2 times the
 * machine epsilon from 1.
 */
TEST_F(BasisSplineTest, BasisSplinePartitionOfUnityTest)
{
    // set maximum degree:
    int max_degree = 5;

    // create basis spline functor:
    BasisSpline B;

    // loop over degrees:
    for(int degree = 0; degree <= max_degree; degree++)
    {
        // get size of basis:
        int nBasis = knotVector_.size() + degree - 1;

        // append and prepend the apprpropriate number of knots:
        std::vector<real> knots;
        for(int i = 0; i < degree; i++)
        {
            knots.push_back(knotVector_.front());
        }
        for(unsigned int i = 0; i < knotVector_.size(); i++)
        {
            knots.push_back(knotVector_[i]);
        }
        for(int i = 0; i < degree; i++)
        {
            knots.push_back(knotVector_.back());
        }
     
        // loop over evaluation points:
        for(unsigned int i = 0; i < evalPoints_.size(); i++)
        {
            // initialise sum as zero:
            real unity = 0.0;

            // loop over basis:
            for(int j = 0; j < nBasis; j++)
            {
                // add value to sum:
                unity += B(knots, degree, j, evalPoints_[i]);
            }

            // assert partition of unity property:
            ASSERT_NEAR(1.0, unity, 2*std::numeric_limits<real>::epsilon());
        }
    }
}


/*!
 * Tests that the BasisSpline functor gives correct values for quadratic 
 * splines. Reference values have been taken from the R software package and 
 * are hardcoded into this test case. Threshold for floating point comparison
 * is taken to be the machine epsilon.
 */
TEST_F(BasisSplineTest, BasisSplineQuadraticTest)
{
    // test third degree / cubic splines:   
    int degree = 2;

    // reference values:
    std::vector<real> refValQuadratic = {1.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 
                                         0.32653060, 0.51275510, 0.16071430, 0.00000000, 0.00000000, 0.00000000,
                                         0.00000000, 0.00000000, 0.50000000, 0.50000000, 0.00000000, 0.00000000, 
                                         0.00000000, 0.00000000, 0.00000000, 0.87500000, 0.12500000, 0.00000000,
                                         0.02040816, 0.33673469, 0.64285714, 0.00000000, 0.00000000, 0.00000000, 
                                         0.00000000, 0.00000000, 0.00000000, 0.47759225, 0.45418029, 0.06822746,
                                         0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.00000000};

    // create BasisSpline functor:
    BasisSpline B;
 
    // append and prepend the apprpropriate number of knots:
    std::vector<real> knots;
    for(int i = 0; i < degree; i++)
    {
        knots.push_back(knotVector_.front());
    }
    for(unsigned int i = 0; i < knotVector_.size(); i++)
    {
        knots.push_back(knotVector_[i]);
    }
    for(int i = 0; i < degree; i++)
    {
        knots.push_back(knotVector_.back());
    }

    // size of basis:
    int nBasis = knotVector_.size() + degree - 1;

    // loop over evalution points:    
    for(unsigned int j = 0; j < evalPoints_.size(); j++)
    {
        // loop ober knot intervals:
        for(int i = 0; i < nBasis; i++)
        {
            real b = B(knots, degree, i, evalPoints_[j]); 
            ASSERT_NEAR(refValQuadratic[j*nBasis + i], b, std::numeric_limits<real>::epsilon());
        }
    }
}


/*!
 * Tests that the BasisSpline functor gives correct values for cubic splines.
 * Reference values have been taken from the R software package and are hard-
 * coded into this test case. Threshold for floating point comparison is taken
 * to be the machine epsilon here.
 */
TEST_F(BasisSplineTest, BasisSplineCubicTest)
{
    // test third degree / cubic splines:   
    int degree = 3;

    // reference values:
    std::vector<real> refValCubic = {1.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
                                     0.18658892, 0.46041363, 0.29942602, 0.05357143, 0.00000000, 0.00000000, 0.00000000,
                                     0.00000000, 0.00000000, 0.05555556, 0.88888889, 0.05555556, 0.00000000, 0.00000000,
                                     0.00000000, 0.00000000, 0.00000000, 0.68055560, 0.30381940, 0.01562500, 0.00000000,
                                     0.00291545, 0.10167639, 0.46683674, 0.42857143, 0.00000000, 0.00000000, 0.00000000,
                                     0.00000000, 0.00000000, 0.00000000, 0.27443368, 0.49676188, 0.21098317, 0.01782128,
                                     0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.00000000};

    // create BasisSpline functor:
    BasisSpline B;

    // append and prepend the apprpropriate number of knots:
    std::vector<real> knots;
    for(int i = 0; i < degree; i++)
    {
        knots.push_back(knotVector_.front());
    }
    for(unsigned int i = 0; i < knotVector_.size(); i++)
    {
        knots.push_back(knotVector_[i]);
    }
    for(int i = 0; i < degree; i++)
    {
        knots.push_back(knotVector_.back());
    }
    
    // size of basis:
    int nBasis = knotVector_.size() + degree - 1;

    // loop over evalution points:    
    for(unsigned int j = 0; j < evalPoints_.size(); j++)
    {
        // loop ober knot intervals:
        for(int i = 0; i < nBasis; i++)
        {      
            real b = B(knots, degree, i, evalPoints_[j]); 
            ASSERT_NEAR(refValCubic[j*nBasis + i], b, std::numeric_limits<real>::epsilon());
        }
    }
}


/*!
 * Tests that the BasisSplineDerivative functor gives correct values for the 
 * first derivative of a cubic spline. Reference values have been taken from 
 * the deriv() function of the splines2 package of the R software. Threshold
 * for floating point comparison is taken to be the machine epsilon here.
 */
TEST_F(BasisSplineTest, BasisSplineFirstDerivativeTest)
{
    // test first derivative here:
    unsigned int derivOrder = 1;

    // test with cubic splines:
    int degree = 3;

    // reference values:
    std::vector<real> refValCubic = {-0.85714286,  0.8571429,  0.0000000,  0.0000000,  0.00000000,  0.0000000, 0.00000000,
                                     -0.27988338, -0.1046829,  0.2774235,  0.1071429,  0.00000000,  0.0000000, 0.00000000,
                                      0.00000000,  0.0000000, -0.3333333,  0.0000000,  0.33333333,  0.0000000, 0.00000000,
                                      0.00000000,  0.0000000,  0.0000000, -0.5833333,  0.48958333,  0.0937500, 0.00000000,
                                     -0.01749271, -0.2350583, -0.1760204,  0.4285714,  0.00000000,  0.0000000, 0.00000000,
                                      0.00000000,  0.0000000,  0.0000000, -0.3183948, -0.02224038,  0.2821545, 0.05848068,
                                      0.00000000,  0.0000000,  0.0000000,  0.0000000,  0.00000000, -0.8571429, 0.85714286};

    // create basis spline derivative functor:
    BasisSplineDerivative dB;

    // append and prepend the appropriate number of knots:
    std::vector<real> knots;
    for(int i = 0; i < degree; i++)
    {
        knots.push_back(knotVector_.front());
    }
    for(unsigned int i = 0; i < knotVector_.size(); i++)
    {
        knots.push_back(knotVector_[i]);
    }
    for(int i = 0; i < degree; i++)
    {
        knots.push_back(knotVector_.back());
    }

    // size of basis:
    int nBasis = knotVector_.size() + degree - 1;

    // loop over evaluation points:    
    for(unsigned int j = 0; j < evalPoints_.size(); j++)
    {
        // loop over knot intervals:
        for(int i = 0; i < nBasis; i++)
        {      
            real d = dB(knots, degree, i, evalPoints_[j], derivOrder);
            ASSERT_NEAR(refValCubic[j*nBasis + i], d, std::numeric_limits<real>::epsilon());
        }
    }
}


/*!
 * Tests that the BasisSplineDerivative functor gives correct values for the 
 * second derivative of a cubic spline. Reference values have been taken from 
 * the deriv() function of the splines2 package of the R software. Threshold
 * for floating point comparison is taken to be the machine epsilon here.
 */
TEST_F(BasisSplineTest, BasisSplineSecondDerivativeTest)
{
    // test first derivative here:
    unsigned int derivOrder = 2;

    // test with cubic splines:
    int degree = 3;

    // reference values:
    std::vector<real> refValCubic = {0.48979592, -0.9183673,  0.42857143,  0.0000000,  0.0000000,  0.00000000, 0.0000000,
                                     0.27988338, -0.3640671, -0.05867347,  0.1428571,  0.0000000,  0.00000000, 0.0000000,
                                     0.00000000,  0.0000000,  1.33333333, -2.6666667,  1.3333333,  0.00000000, 0.0000000,
                                     0.00000000,  0.0000000,  0.00000000,  0.3333333, -0.7083333,  0.37500000, 0.0000000,
                                     0.06997085,  0.1902332, -0.54591837,  0.2857143,  0.0000000,  0.00000000, 0.0000000,
                                     0.00000000,  0.0000000,  0.00000000,  0.2462654, -0.4113694,  0.03716744, 0.1279366,
                                     0.00000000,  0.0000000,  0.00000000,  0.0000000,  0.4285714, -0.91836735, 0.4897959};

    // create basis spline derivative functor:
    BasisSplineDerivative dB;

    // append and prepend the apprpropriate number of knots:
    std::vector<real> knots;
    for(int i = 0; i < degree; i++)
    {
        knots.push_back(knotVector_.front());
    }
    for(unsigned int i = 0; i < knotVector_.size(); i++)
    {
        knots.push_back(knotVector_[i]);
    }
    for(int i = 0; i < degree; i++)
    {
        knots.push_back(knotVector_.back());
    }

    // size of basis:
    int nBasis = knotVector_.size() + degree - 1;

    // loop over evalution points:    
    for(unsigned int j = 0; j < evalPoints_.size(); j++)
    {
        // loop ober knot intervals:
        for(int i = 0; i < nBasis; i++)
        {      
            real d = dB(knots, degree, i, evalPoints_[j], derivOrder);
            ASSERT_NEAR(refValCubic[j*nBasis + i], d, std::numeric_limits<real>::epsilon());
        }
    }
}

