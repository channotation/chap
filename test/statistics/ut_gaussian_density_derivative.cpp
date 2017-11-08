#include <limits>
#include <random>

#include <gtest/gtest.h>

#include "statistics/gaussian_density_derivative.hpp"


/*!
 * \brief Test fixture for the GaussianDensityDerivative.
 */
class GaussianDensityDerivativeTest : public ::testing::Test
{
    public:

        /*!
         * Constructor is used to set up a random sample drawn from as 
         * Gaussian distribution.
         */
        GaussianDensityDerivativeTest()
        {
            // parameters of normal distribution:
            real mu = -1.0;
            real sd = 3.0;

            // prepare random distribution:
            std::default_random_engine generator;
            std::normal_distribution<real> distribution(mu, sd);

            // create a random sample:
            size_t numSamples = 1e2;
            for(size_t i = 0; i < numSamples; i++)
            {
                testData_.push_back( distribution(generator) );
            }
        };

    protected:

        std::vector<real> testData_;
};


/*!
 *
 */
TEST_F(GaussianDensityDerivativeTest, GaussianDensityDerivativeCoefATest)
{
    real tolerance = std::numeric_limits<real>::epsilon();

    // prepare density derivative estimator:
    GaussianDensityDerivative gdd;
    gdd.setDerivOrder(2);
    gdd.setBandWidth(0.1);
    gdd.setErrorBound(0.001);

    // calculate coefficients using reference and optimised implementation:
    std::vector<real> coefA = gdd.setupCoefA();
//    std::vector<real> coefARef = gdd.setupCoefARef();
/*
    // both coefficient vectors should have same size:
    ASSERT_EQ(coefA.size(), coefARef.size());

    // check that both functions return same values:
    for(int i = 0; i < coefARef.size(); i++)
    {
        ASSERT_NEAR(coefARef[i], coefA[i], tolerance);   
    }*/
}


/*!
 *
 */
TEST_F(GaussianDensityDerivativeTest, GaussianDensityDerivativeConsistencyTest)
{
    std::vector<real> sample = testData_;
/*

    std::vector<real> eval = {0.0};

    real eps = 0.01;

    
    GaussianDensityDerivative gdd;
    gdd.setDerivOrder(2);
    gdd.setBandWidth(0.1);
    gdd.setErrorBound(eps);





    // estimate derivative via direct loop and via approximate method:
    std::vector<real> derivDirect = gdd.estimateDirect(sample, eval);
    std::vector<real> derivApprox = gdd.estimateApprox(sample, eval);
    
    // loop over all eval points and assert equality of estimation methods:
    for(int i = 0; i < eval.size(); i++)
    {
        std::cout<<"derivDirect = "<<derivDirect[i]<<"  "
                 <<"derivApprox = "<<derivApprox[i]<<"  "
                 <<std::endl;
        ASSERT_NEAR(derivDirect[i], derivApprox[i], eps);
    }*/
}

