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


#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "geometry/bspline_basis_set.hpp"


/*!
 * \brief Test fixture for BSplineBasisSet functor.
 */
class BSplineBasisSetTest : public ::testing::Test
{
		protected:

        // knot vector:
        std::vector<real> uniqueKnots_ = {-4, -0.5, 0.0, 0.5, 4};
         
        // evaluation points:
        std::vector<real> evalPoints_ = {-4.0, -2.5, 0.0, 0.5, -1.0, std::sqrt(2.0), 4.0};

        // true knot span index:
        std::vector<size_t> knotSpanIdx_ = {0, 0, 2, 3, 0, 3, 3};

        // create degree-appropriate knot vector from unique knots:
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
 * Test that BSplineBasisSet finds the correct knot span for a number of
 * evaluation points.
 */
TEST_F(BSplineBasisSetTest, BSplineBasisSetKnotSpanTest)
{
    // create functor:
    BSplineBasisSet B;

    // set degree and adapt knotVector approproately:
    unsigned int degree = 3;
    std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints_.size(); i++)
    {
        // asser correct knot span index:
        size_t idx = B.findKnotSpan(evalPoints_[i], knots, degree);
        ASSERT_EQ(knotSpanIdx_[i] + degree, idx);
    }
}


/*!
 * Checks that the sum over the entire basis is equal to one for constant 
 * degree. This is done for spline of up to fifths degree.
 */
TEST_F(BSplineBasisSetTest, BSplineBasisSetParitionOfUnityTest)
{ 
    // create spline basis functor:
    BSplineBasisSet B;

    // loop over various degrees:
    unsigned int maxDegree = 5;
    for(unsigned int degree = 0; degree <= maxDegree; degree++)
    {
        // prepare knots for this degree:
        std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

        // loop over evaluation points:
        for(size_t i = 0; i < evalPoints_.size(); i++)
        {
            // evaluate basis:
            SparseBasis basis = B(evalPoints_[i], knots, degree);

            // compute sum over basis functions:
            real unity = 0.0;
            for(auto b = basis.begin(); b != basis.end(); b++)
            {
                unity += b -> second;
            }

            // basis should sum to one:
            ASSERT_NEAR(1.0, unity, 10*std::numeric_limits<real>::epsilon());
        }
    }
}


/*!
 * Checks that the BSplineBasisSet functor returns correct values for quadratic
 * splines. This is done by asserting agreement with reference values computed
 * with the R software, which are hardcoded into this test case. The threshold
 * for floating point comparison is taken to be the machine precision here.
 */
TEST_F(BSplineBasisSetTest, BSplineBasisSetQuadraticTest)
{
    // test second degree / quadratic splines:   
    unsigned int degree = 2;

    // reference values:
    int nBasis = 6;
    std::vector<real> refVal = {
            1.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 
            0.32653060, 0.51275510, 0.16071430, 0.00000000, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.50000000, 0.50000000, 0.00000000, 0.00000000, 
            0.00000000, 0.00000000, 0.00000000, 0.87500000, 0.12500000, 0.00000000,
            0.02040816, 0.33673469, 0.64285714, 0.00000000, 0.00000000, 0.00000000, 
            0.00000000, 0.00000000, 0.00000000, 0.47759225, 0.45418029, 0.06822746,
            0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.00000000};

    // create basis set functor:
    BSplineBasisSet B;
 
    // append and prepend the apprpropriate number of knots:
    std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints_.size(); i++)
    {
        // evaluate basis set:
        SparseBasis basis = B(evalPoints_[i], knots, degree);

        // loop over basis:
        for(auto b : basis)
        {
            // check agreement with reference values:
            ASSERT_NEAR(
                    refVal[i*nBasis + b.first], 
                    b.second,
                    std::numeric_limits<real>::epsilon());
        }
    }
}


/*!
 * Checks that the BSplineBasisSet functor returns correct values for cubic 
 * splines. This is done by asserting agreement with reference values computed
 * with the R software, which are hardcoded into this test case. The threshold
 * for floating point comparison is taken to be the machine precision here.
 */
TEST_F(BSplineBasisSetTest, BSplineBasisSetCubicTest)
{
    // test third degree / cubic splines:   
    unsigned int degree = 3;

    // reference values:
    int nBasis = 7;
    std::vector<real> refVal = {
            1.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
            0.18658892, 0.46041363, 0.29942602, 0.05357143, 0.00000000, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.05555556, 0.88888889, 0.05555556, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.00000000, 0.68055560, 0.30381940, 0.01562500, 0.00000000,
            0.00291545, 0.10167639, 0.46683674, 0.42857143, 0.00000000, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.00000000, 0.27443368, 0.49676188, 0.21098317, 0.01782128,
            0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.00000000};

    // create basis set functor:
    BSplineBasisSet B;
 
    // append and prepend the apprpropriate number of knots:
    std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints_.size(); i++)
    {
        // evaluate basis set:
        SparseBasis basis = B(evalPoints_[i], knots, degree);

        // loop over basis:
        for(auto b : basis)
        {
            // check agreement with reference values:
            ASSERT_NEAR(
                    refVal[i*nBasis + b.first], 
                    b.second,
                    std::numeric_limits<real>::epsilon());
        }
    }
}


/*!
 * Checks that the functor returns the correct values of the B-spline basis if
 * evaluated through the derivative interface.
 */
TEST_F(BSplineBasisSetTest, BSplineBasisSetZerothDerivativeTest)
{
    // specify degree and derivative order:
    int deriv = 0;
    int degree = 3;

    // reference values:
    int nBasis = 7;
    std::vector<real> refVal = {
            1.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
            0.18658892, 0.46041363, 0.29942602, 0.05357143, 0.00000000, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.05555556, 0.88888889, 0.05555556, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.00000000, 0.68055560, 0.30381940, 0.01562500, 0.00000000,
            0.00291545, 0.10167639, 0.46683674, 0.42857143, 0.00000000, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.00000000, 0.27443368, 0.49676188, 0.21098317, 0.01782128,
            0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.00000000};

    // create basis set functor:
    BSplineBasisSet B;

    // create appropriate knot vector:
    std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints_.size(); i++)
    {
        // evaluate basis (derivatives) at this point:
        SparseBasis basis = B(
                evalPoints_[i],
                knots,
                degree,
                deriv);

        // for zeroth order derivative will perform partition of unity test:
        real unity = 0.0;

        // loop over basis (derivatives):
        for(auto b : basis)
        {
            // check agreement with reference values:            
            ASSERT_NEAR(
                    refVal[i*nBasis + b.first],
                    b.second,
                    std::numeric_limits<real>::epsilon());

            // increment sum over basis elements:
            unity += b.second;
        }

        // check partition of unity property:
        ASSERT_NEAR(1.0, unity, std::numeric_limits<real>::epsilon());
    }
}


/*!
 * Checks that the BSplineBasisSet functor returns the correct values of the 
 * first derivative for a given set of test points.
 */
TEST_F(BSplineBasisSetTest, BSplineBasisSetFirstDerivativeTest)
{
    // specify degree and derivative order:
    int deriv = 1;
    int degree = 3;

    // reference values:
    int nBasis = 7;
    std::vector<real> refVal = {
            -0.85714286,  0.8571429,  0.0000000,  0.0000000,  0.00000000,  0.0000000, 0.00000000,
            -0.27988338, -0.1046829,  0.2774235,  0.1071429,  0.00000000,  0.0000000, 0.00000000,
             0.00000000,  0.0000000, -0.3333333,  0.0000000,  0.33333333,  0.0000000, 0.00000000,
             0.00000000,  0.0000000,  0.0000000, -0.5833333,  0.48958333,  0.0937500, 0.00000000,
            -0.01749271, -0.2350583, -0.1760204,  0.4285714,  0.00000000,  0.0000000, 0.00000000,
             0.00000000,  0.0000000,  0.0000000, -0.3183948, -0.02224038,  0.2821545, 0.05848068,
             0.00000000,  0.0000000,  0.0000000,  0.0000000,  0.00000000, -0.8571429, 0.85714286};

    // create basis set functor:
    BSplineBasisSet B;

    // create appropriate knot vector:
    std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints_.size(); i++)
    {
        // evaluate basis (derivatives) at this point:
        SparseBasis basis = B(
                evalPoints_[i],
                knots,
                degree,
                deriv);

        // loop over basis (derivatives):
        for(auto b : basis)
        {
            // check agreement with reference values:
            ASSERT_NEAR(
                    refVal[i*nBasis + b.first],
                    b.second,
                    std::numeric_limits<real>::epsilon());
        }
    }
}


/*!
 * Checks that the BSplineBasisSet functor returns the correct values for the
 * second derivatives of the basis functions for a given set of test points.
 */
TEST_F(BSplineBasisSetTest, BSplineBasisSetSecondDerivativeTest)
{
    // specify degree and derivative order:
    int deriv = 2;
    int degree = 3;

    // reference values:
    int nBasis = 7;
    std::vector<real> refVal = {
            0.48979592, -0.9183673,  0.42857143,  0.0000000,  0.0000000,  0.00000000, 0.0000000,
            0.27988338, -0.3640671, -0.05867347,  0.1428571,  0.0000000,  0.00000000, 0.0000000,
            0.00000000,  0.0000000,  1.33333333, -2.6666667,  1.3333333,  0.00000000, 0.0000000,
            0.00000000,  0.0000000,  0.00000000,  0.3333333, -0.7083333,  0.37500000, 0.0000000,
            0.06997085,  0.1902332, -0.54591837,  0.2857143,  0.0000000,  0.00000000, 0.0000000,
            0.00000000,  0.0000000,  0.00000000,  0.2462654, -0.4113694,  0.03716744, 0.1279366,
            0.00000000,  0.0000000,  0.00000000,  0.0000000,  0.4285714, -0.91836735, 0.4897959};

    // create basis set functor:
    BSplineBasisSet B;

    // create appropriate knot vector:
    std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints_.size(); i++)
    {
        // evaluate basis (derivatives) at this point:
        SparseBasis basis = B(
                evalPoints_[i],
                knots,
                degree,
                deriv);

        // loop over basis (derivatives):
        for(auto b : basis)
        {
            // check agreement with reference values:
            ASSERT_NEAR(
                    refVal[i*nBasis + b.first],
                    b.second,
                    std::numeric_limits<real>::epsilon());
        }
    }
}

