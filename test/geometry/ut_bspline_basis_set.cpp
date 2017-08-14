#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>

#include <gtest/gtest.h>

#include "geometry/bspline_basis_set.hpp"


/*
 *
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
            for(int i = 0; i < degree; i++)
            {
                knots.push_back(uniqueKnots.front());
            }
            for(unsigned int i = 0; i < uniqueKnots.size(); i++)
            {
                knots.push_back(uniqueKnots[i]);
            }
            for(int i = 0; i < degree; i++)
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
    unsigned int degree = 2;
    std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints_.size(); i++)
    {
        // asser correct knot span index:
        size_t idx = B.findKnotSpan(evalPoints_[i], knots, degree);
        ASSERT_EQ(knotSpanIdx_[i] + degree, idx);
    }
}


/*
 *
 */
TEST_F(BSplineBasisSetTest, BSplineBasisSetParitionOfUnityTest)
{
    
    // create spline basis functor:
    BSplineBasisSet B;

    // evaluate only spline itself, not derivative:
    int deriv = 0;

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
            std::vector<real> basis = B(evalPoints_[i], knots, degree, deriv);

            // compute sum over basis functions:
            real unity = 0.0;
            for(size_t j = 0; j < basis.size(); j++)
            {
    /*            std::cout<<"j = "<<j<<"  "
                         <<"basis[j] = "<<basis[j]<<"  "
                         <<std::endl;*/

                unity += basis[j];
            }

            // basis should sum to one:
            ASSERT_NEAR(1.0, unity, std::numeric_limits<real>::epsilon());
        }
    }
}



































