#include <vector>
#include <cmath>

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
        std::vector<real> evalPoints_ = {-4.0, -2.5, 0.0, 0.5, -1.0, std::sqrt(2.0)};
};


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
                                     0.00000000, 0.00000000, 0.00000000, 0.27443368, 0.49676188, 0.21098317, 0.01782128};

    // create BasisSpline functor:
    BasisSpline bs;
   
    // size of basis:
    int nBasis = knotVector_.size() - degree - 1;

    // loop over evalution points:    
    for(int j = 0; j < evalPoints_.size(); j++)
    {
        // loop ober knot intervals:
        for(int i = 0; i < nBasis; i++)
        {      
            real b = bs(knotVector_, degree, i, evalPoints_[j]); 

            /*
            std::cout<<"B_"<<i<<"_"<<degree<<"("<<evalPoints_[j]<<")"<<" = "
                     <<b<<std::endl; 
            std::cout<<"refVal = "<<refValCubic[j*nBasis + i]<<std::endl;
            */

            ASSERT_NEAR(refValCubic[j*nBasis + i], b, 1e-7);
        }
    }
}

