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
 * Checks that the sum over all basis functions for the entire basis set equals
 * one for fixed degree. This is assessed for splines up to degree five.
 */
TEST_F(BSplineBasisElementTest, BSplineBasisElementPartitionOfUnityTest)
{
    // create functor:
    BSplineBasisElement B;

    int maxDegree = 5;

    // loop over degrees:
    for(int degree = 0; degree <= maxDegree; degree++)
    {
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
                unity += basisElement;
            }

            // assert partition of unity property:
            ASSERT_NEAR(1.0, unity, std::numeric_limits<real>::epsilon());
        }
    }
}


/*!
 * Checks that the BSplineBasisElement functor gives the correct values for a
 * quadratic spline. This is assessed by comparison to reference values 
 * computes with the R software. The floating point comparison threshold is
 * the machine epsilon.
 */
TEST_F(BSplineBasisElementTest, BSplineBasisElementQuadraticTest)
{
    // test second degree / quadratic splines:   
    int degree = 2;

    // reference values:
    std::vector<real> refVal = {
            1.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 
            0.32653060, 0.51275510, 0.16071430, 0.00000000, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.50000000, 0.50000000, 0.00000000, 0.00000000, 
            0.00000000, 0.00000000, 0.00000000, 0.87500000, 0.12500000, 0.00000000,
            0.02040816, 0.33673469, 0.64285714, 0.00000000, 0.00000000, 0.00000000, 
            0.00000000, 0.00000000, 0.00000000, 0.47759225, 0.45418029, 0.06822746,
            0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.00000000};

    // create BasisSpline functor:
    BSplineBasisElement B;

    // set up knot vector:
    std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

    // size of complete basis:
    size_t nBasis = knots.size() - degree - 1;

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints_.size(); i++)
    {
        // loop over basis elements:
        for(size_t j = 0; j < nBasis; j++)
        {
            real basisElement = B(evalPoints_[i], j, knots, degree);
            ASSERT_NEAR(
                    refVal[i*nBasis + j], 
                    basisElement, 
                    std::numeric_limits<real>::epsilon());
        }
    }
}


/*!
 * Checks that the BSplineBasisElement functor gives the correct values for a
 * cubic spline. This is assessed by comparison to reference values 
 * computes with the R software. The floating point comparison threshold is
 * the machine epsilon.
 */
TEST_F(BSplineBasisElementTest, BSplineBasisElementCubicTest)
{
    // test third degree / cubic splines:   
    int degree = 3;

    // reference values:
    std::vector<real> refVal = {
            1.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
            0.18658892, 0.46041363, 0.29942602, 0.05357143, 0.00000000, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.05555556, 0.88888889, 0.05555556, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.00000000, 0.68055560, 0.30381940, 0.01562500, 0.00000000,
            0.00291545, 0.10167639, 0.46683674, 0.42857143, 0.00000000, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.00000000, 0.27443368, 0.49676188, 0.21098317, 0.01782128,
            0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.00000000};

    // create BasisSpline functor:
    BSplineBasisElement B;

    // set up knot vector:
    std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

    // size of complete basis:
    size_t nBasis = knots.size() - degree - 1;

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints_.size(); i++)
    {
        // loop over basis elements:
        for(size_t j = 0; j < nBasis; j++)
        {
            real basisElement = B(evalPoints_[i], j, knots, degree);
            ASSERT_NEAR(
                    refVal[i*nBasis + j], 
                    basisElement, 
                    std::numeric_limits<real>::epsilon());
        }
    }
}


/*
 *
 */
TEST_F(BSplineBasisElementTest, BSplineBasisElementZerothDerivativeTest)
{
    // specify degree and derivative order:
    int deriv = 0;
    int degree = 3;

    // reference values:
    std::vector<real> refVal = {
            1.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
            0.18658892, 0.46041363, 0.29942602, 0.05357143, 0.00000000, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.05555556, 0.88888889, 0.05555556, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.00000000, 0.68055560, 0.30381940, 0.01562500, 0.00000000,
            0.00291545, 0.10167639, 0.46683674, 0.42857143, 0.00000000, 0.00000000, 0.00000000,
            0.00000000, 0.00000000, 0.00000000, 0.27443368, 0.49676188, 0.21098317, 0.01782128,
            0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 1.00000000};

    // create basis set functor:
    BSplineBasisElement B;

    // create appropriate knot vector:
    std::vector<real> knots = prepareKnotVector(uniqueKnots_, degree);

    // number of basis functions:
    int nBasis = knots.size() - degree - 1;

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints_.size(); i++)
    {
        // for zeroth order derivative will perform partition of unity test:
        real unity = 0.0;

        // loop over basis (derivatives):
        for(int j = 0; j < nBasis; j++)
        {
            // evaluate this basis function:
            real basisElement = B(evalPoints_[i], knots, degree, j, deriv);

            std::cout<<"eval = "<<evalPoints_[i]<<"  "
                     <<"j = "<<j<<"  "
                     <<"B = "<<basisElement<<"  "
                     <<std::endl;


            // check agreement with reference values:            
            ASSERT_NEAR(
                    refVal[i*nBasis + j],
                    basisElement,
                    std::numeric_limits<real>::epsilon());

            // increment sum over basis elements:
            unity += basisElement;
        }

        std::cout<<std::endl;

        // check partition of unity property:
        ASSERT_NEAR(1.0, unity, std::numeric_limits<real>::epsilon());
    }
}







