#include <algorithm>
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
            size_t numSamples = 1e1;
            for(size_t i = 0; i < numSamples; i++)
            {
                testData_.push_back( distribution(generator) );
            }
        };

    protected:

        std::vector<real> testData_;
};

/*!
 * Checks that space partitioning produces correct centres and that each data
 * point is associated with the correct centre.
 */
TEST_F(GaussianDensityDerivativeTest, 
       GaussianDensityDerivativeSpacePartitioningTest)
{
    real tolerance = std::numeric_limits<real>::epsilon();

    // run this test for a variaty of bandwidths:
    std::vector<real> bandwidth = {2.0, 1.0, 0.5, 0.3, 0.1, 0.0001};
    for(auto bw : bandwidth)
    {
        // prepare density derivative estimator:
        GaussianDensityDerivative gdd;
        gdd.setDerivOrder(2);
        gdd.setBandWidth(bw);
        gdd.setErrorBound(0.001);

        // obtain cluster centres:
        std::vector<real> centres = gdd.setupClusterCentres();

        // largest and smallest cluster centre:
        real cMin = *std::min_element(centres.begin(), centres.end());
        real cMax = *std::max_element(centres.begin(), centres.end());

        // check that all centres fall into unit interval:
        ASSERT_GE(cMax, cMin);
        ASSERT_GE(cMin, 0.0);
        ASSERT_GE(1.0, cMax);

        // check that centre distance is bounded by half bandwidth: 
        for(size_t i = 0; i < centres.size() - 1; i++)
        {
            ASSERT_LE(centres[i+1] - centres[i], bw/2.0 + tolerance);
        }

        // sample data in interval [0,1]:
        std::vector<real> sample = {0.2, 0.1, 0.33, 0.4, 0.1, 0.9999, 0.7};

        // find interval centre for each data point:
        std::vector<unsigned int> idx = gdd.setupClusterIndices(sample);

        // for each data point check that correct centre was found: 
        for(size_t i = 0; i < sample.size(); i++)
        {
            // find distance from closest centre on record:
            real d = sample.at(i) - centres.at(idx.at(i));
            
            // calculate distance from sample to all centres:
            std::vector<real> dist;
            dist.reserve(centres.size());
            for(auto c : centres)
            {
                dist.push_back(std::abs(sample[i] - c));
            }

            // find minimal distance:
            auto minIter = min_element(dist.begin(), dist.end());
            real minDist = *minIter;
            real minIdx = minIter - dist.begin();

            // check that centre on record is actually closest:
            ASSERT_LE(d, minDist);
        }
    }
}


/*!
 * Checks that the a-coefficients are correctly computed for two different
 * derivative orders by comparison to manually calculated coefficients.
 */
TEST_F(GaussianDensityDerivativeTest, GaussianDensityDerivativeCoefATest)
{
    real tolerance = std::numeric_limits<real>::epsilon();

    // prepare density derivative estimator:
    GaussianDensityDerivative gdd;
    gdd.setDerivOrder(2);
    gdd.setBandWidth(0.1);
    gdd.setErrorBound(0.001);

    // manually computed coefficients for this parameter vector:
    std::vector<real> coefATrue = {1, -2, 1, -1};

    // calculate coefficients using reference and optimised implementation:
    std::vector<real> coefA = gdd.setupCoefA();

    // check correct number of coefficients:
    ASSERT_EQ(coefATrue.size(), coefA.size());

    // check correct value of coefficients:
    for(size_t i = 0; i < coefATrue.size(); i++)
    {
        ASSERT_NEAR(coefATrue[i], coefA[i], tolerance);
    }

    // set a different derivative order:
    gdd.setDerivOrder(5);

    // manually computed coefficients for this parameter vector:
    coefATrue = {  1, -5,  10, -10,  5,  -1, 
                 -10, 30, -30,  10, 15, -15};

    // calculate coefficients using reference and optimised implementation:
    coefA = gdd.setupCoefA();

    // check correct number of coefficients:
    ASSERT_EQ(coefATrue.size(), coefA.size());

    // check correct value of coefficients:
    for(size_t i = 0; i < coefATrue.size(); i++)
    {
        ASSERT_NEAR(coefATrue[i], coefA[i], tolerance);
    }
}


/*!
 */
TEST_F(GaussianDensityDerivativeTest, GaussianDensityDerivativeCoefBTest)
{
    real tolerance = std::numeric_limits<real>::epsilon();

    // prepare density derivative estimator:
    unsigned int deriv = 2;
    GaussianDensityDerivative gdd;
    gdd.setDerivOrder(deriv);
    gdd.setBandWidth(0.1);
    gdd.setErrorBound(0.001);

    // sample data in interval [0,1]:
    std::vector<real> sample = {0.0, 0.33, 0.5, 1.0};


    // prepare spatial partitioning:
    gdd.centres_ = gdd.setupClusterCentres();
    gdd.idx_ = gdd.setupClusterIndices(sample);

    // calculate parameters need for coefficient calculation:
    gdd.q_ = gdd.setupCoefQ(sample.size());
    gdd.epsPrime_ = gdd.setupScaledTolerance(sample.size());
    gdd.rc_ = gdd.setupCutoffRadius();
    gdd.trunc_ = gdd.setupTruncationNumber();

    
    std::vector<real> coefBTrue(gdd.trunc_*(deriv+1)*gdd.centres_.size(), 0.0);


    // calculate coefficients using reference and optimised implementation:
    std::vector<real> coefB = gdd.setupCoefB(sample);
    std::vector<real> b = gdd.compute_B(sample);

    // check right number of coefficients:
    ASSERT_EQ(coefBTrue.size(), coefB.size());


    std::cout<<"b.size = "<<b.size()<<std::endl;
    std::cout<<"coefB.size = "<<coefB.size()<<std::endl;

    // assert that correct B coefficients have been computed:
    for(size_t i = 0; i < b.size(); i++)
    {
        std::cout<<"B_true = "<<coefBTrue[i]<<"  "
                 <<"coefB = "<<coefB[i]<<"  "
                 <<"B = "<<b[i]<<"  "
                 <<std::endl;

//        ASSERT_NEAR(coefBTrue[i], coefB[i], tolerance);
    }
}


/*!
 *
 */
TEST_F(GaussianDensityDerivativeTest, GaussianDensityDerivativeConsistencyTest)
{

    std::vector<real> sample = {0.0, 1.0};

    std::vector<real> eval = sample;

    real eps = 0.001;

    
    GaussianDensityDerivative gdd;
    gdd.setDerivOrder(2);
    gdd.setBandWidth(1.0);
    gdd.setErrorBound(eps);





    // estimate derivative via direct loop and via approximate method:
    std::vector<real> derivDirect = gdd.estimateDirect(sample, eval);
    std::vector<real> derivApprox = gdd.estimateApprox(sample, eval);
    
    // loop over all eval points and assert equality of estimation methods:
    for(int i = 0; i < eval.size(); i++)
    {
//        std::cout<<"derivDirect = "<<derivDirect[i]<<"  "
//                 <<"derivApprox = "<<derivApprox[i]<<"  "
//                 <<std::endl;
//        ASSERT_NEAR(derivDirect[i], derivApprox[i], eps);
    }
}

