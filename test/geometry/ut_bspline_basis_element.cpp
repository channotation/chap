#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>

#include <gtest/gtest.h>

#include "geometry/bspline_basis_element.hpp"


/*
 *
 */
class BSplineBasisElementTest : public ::testing::Test
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
 *
 */
TEST_F(BSplineBasisElementTest, BSplineBasisElementPartitionOfUnityTest)
{
    // create functor:
    BSplineBasisElement B;

    int degree = 3;

    // set up knot vector:
    std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

    // size of complete basis vector:
    int nBasis = knots.size() - degree - 1;

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints_.size(); i++)
    {
        real unity = 0.0;

        for(size_t j = 0; j < nBasis; j++)
        {
            real basisElement = B(evalPoints_[i], j, knots, degree);
            std::cout<<"eval = "<<evalPoints_[i]<<"  "
                     <<"j = "<<j<<"  "
                     <<"B = "<<basisElement<<"  "
                     <<std::endl;
            unity += basisElement;
        }

        // assert partition of unity property:
        ASSERT_NEAR(1.0, unity, std::numeric_limits<real>::epsilon());
    }
}














