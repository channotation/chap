#include <vector>
#include <cmath>
#include <iostream>
#include <limits>

#include <gtest/gtest.h>

#include "geometry/basis_spline.hpp"


/*
 * Test fixture for testing the BasisSpline class. Reference values have been
 * computed with the bs() method of the spline library in R.
 */
class BasisSplineTest : public ::testing::Test
{
	protected:

        // knot vector:
        std::vector<real> knotVector_ = {-4, -4, -4, -4, -0.5, 0.0, 0.5, 4, 4, 4, 4};
         
        // evaluation points:
        std::vector<real> evalPoints_ = {-4.0, -2.5, 0.0, 0.5, -1.0, std::sqrt(2.0), 4.0};
};


/*
 * Tests that the basis splines over a knot vector form a partition of unity,
 * i.e. that sum_i=1^n B_i,k(x) = 1.
 */
TEST_F(BasisSplineTest, BasisSplinePartitionOfUnityTest)
{
    // set cubic degree:
    int degree = 3;

    // create basis spline functor:
    BasisSpline B;

    // get size of basis:
    int nBasis = knotVector_.size() - degree - 1;

//    evalPoints_ = {4.0};

    // loop over evaluation points:
    for(int i = 0; i < evalPoints_.size(); i++)
    {
        // initialise sum as zero:
        real unity = 0.0;

        // loop over basis:
        for(int j = 0; j < nBasis; j++)
        {/*
            std::cout<<"j = "<<j<<"  "
                     <<"degree = "<<degree_<<"  "
                     <<"unity = "<<unity<<"  "
                     <<"B = "<<B(knotVector_, degree, j, evalPoints_[i])<<std::endl;
     */
            // add value to sum:
            unity += B(knotVector_, degree, j, evalPoints_[i]);
        }

        // assert partition of unity property:
        ASSERT_NEAR(1.0, unity, std::sqrt(std::numeric_limits<real>::epsilon()));
    }
}


/*
 * Tests that the BasisSpline functor gives correct values for quadratic 
 * splines. Reference values have been taken from the R software package and 
 * are hardcoded into this test case.
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
  
    knotVector_ = { -4, -4, -4, -0.5, 0.0, 0.5, 4, 4, 4};

    // size of basis:
    int nBasis = knotVector_.size() - degree - 1;
    std::cout<<"nBasis = "<<nBasis<<std::endl;

    // loop over evalution points:    
    for(int j = 0; j < evalPoints_.size(); j++)
    {
        // loop ober knot intervals:
        for(int i = 0; i < nBasis; i++)
        {      
            real b = B(knotVector_, degree, i, evalPoints_[j]); 

           /* 
            std::cout<<"B_"<<i<<"_"<<degree<<"("<<evalPoints_[j]<<")"<<" = "
                     <<b<<std::endl; 
            std::cout<<"refVal = "<<refValQuadratic[j*nBasis + i]<<std::endl;
            */

            ASSERT_NEAR(refValQuadratic[j*nBasis + i], b, 1e-7);
        }
    }
}


/*
 * Tests that the BasisSpline functor gives correct values for cubic splines.
 * Reference values have been taken from the R software package and are hard-
 * coded into this test case.
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
   
    // size of basis:
    int nBasis = knotVector_.size() - degree - 1;
    std::cout<<"nBasis = "<<nBasis<<std::endl;

    // loop over evalution points:    
    for(int j = 0; j < evalPoints_.size(); j++)
    {
        // loop ober knot intervals:
        for(int i = 0; i < nBasis; i++)
        {      
            real b = B(knotVector_, degree, i, evalPoints_[j]); 

            /*
            std::cout<<"B_"<<i<<"_"<<degree<<"("<<evalPoints_[j]<<")"<<" = "
                     <<b<<std::endl; 
            std::cout<<"refVal = "<<refValCubic[j*nBasis + i]<<std::endl;
            */

            ASSERT_NEAR(refValCubic[j*nBasis + i], b, 1e-7);
        }
    }
}

