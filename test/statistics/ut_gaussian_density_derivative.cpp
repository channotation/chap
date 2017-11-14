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
        // sample data in interval [0,1]:
        std::vector<real> sample = {0.0, 0.1, 0.33, 0.4, 0.1, 0.9999, 0.7, 1.0};
        
        // prepare density derivative estimator:
        GaussianDensityDerivative gdd;
        gdd.setDerivOrder(2);
        gdd.setBandWidth(bw);
        gdd.setErrorBound(0.001);

        gdd.q_ = gdd.setupCoefQ(sample.size());
        gdd.epsPrime_ = gdd.setupScaledTolerance(sample.size());
        gdd.rc_ = gdd.setupCutoffRadius();
        gdd.trunc_ = gdd.setupTruncationNumber();

        // obtain cluster centres:
        std::vector<real> centres = gdd.setupClusterCentres();
        gdd.centres_ = centres;

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
    real tolerance = 1e2*std::numeric_limits<real>::epsilon();

    // prepare density derivative estimator:
    unsigned int deriv = 2;
    GaussianDensityDerivative gdd;
    gdd.setDerivOrder(deriv);
    gdd.setBandWidth(0.1);
    gdd.setErrorBound(0.001);

    // sample data in interval [0,1]:
    std::vector<real> sample = {0.0, 1.0};


    // prepare spatial partitioning:
    gdd.centres_ = gdd.setupClusterCentres();
    gdd.idx_ = gdd.setupClusterIndices(sample);

    // calculate parameters need for coefficient calculation:
    gdd.q_ = gdd.setupCoefQ(sample.size());
    gdd.epsPrime_ = gdd.setupScaledTolerance(sample.size());
    gdd.rc_ = gdd.setupCutoffRadius();
    gdd.trunc_ = gdd.setupTruncationNumber();

    // reference values for coefficients:
    int nCoef = gdd.trunc_*(deriv+1)*gdd.centres_.size(); 
    std::vector<real> coefBTrue(nCoef, 0.0);
    coefBTrue[ 0] = 193.3340454;
    coefBTrue[ 1] = -48.33351135;
    coefBTrue[ 2] =  12.08337784;
    coefBTrue[ 3] = -48.3335113;
    coefBTrue[ 4] =  12.08337784;
    coefBTrue[ 5] =  -3.02084446;
    coefBTrue[ 6] =   6.041688919;
    coefBTrue[ 7] =  -1.51042223;
    coefBTrue[ 8] =   0.3776055574;
    coefBTrue[ 9] =  -0.5034740567;
    coefBTrue[10] =   0.1258685142;
    coefBTrue[11] =  -0.03146712855;
    coefBTrue[12] =   0.03146712855;
    coefBTrue[13] =  -0.007866782136;
    coefBTrue[14] =   0.001966695534;
    coefBTrue[15] =  -0.001573356567;
    coefBTrue[16] =   0.0003933391417;
    coefBTrue[17] =  -9.833478543e-05;
    coefBTrue[18] =   6.555652362e-05;
    coefBTrue[19] =  -1.638913091e-05;
    coefBTrue[20] =   4.097282726e-06;
    coefBTrue[nCoef - 21] = 193.3340302;
    coefBTrue[nCoef - 20] =  48.3335762;
    coefBTrue[nCoef - 19] =  12.08341122;
    coefBTrue[nCoef - 18] =  48.3335762;
    coefBTrue[nCoef - 17] =  12.08341122;
    coefBTrue[nCoef - 16] =   3.020857096;
    coefBTrue[nCoef - 15] =   6.041705608;
    coefBTrue[nCoef - 14] =   1.510428548;
    coefBTrue[nCoef - 13] =   0.3776076734;
    coefBTrue[nCoef - 12] =   0.5034762025;
    coefBTrue[nCoef - 11] =   0.5034762025;
    coefBTrue[nCoef - 10] =   0.03146735206;
    coefBTrue[nCoef -  9] =   0.03146730736;
    coefBTrue[nCoef -  8] =   0.007866838016;
    coefBTrue[nCoef -  7] =   0.001966712065;
    coefBTrue[nCoef -  6] =   0.00157336751;
    coefBTrue[nCoef -  5] =   0.0003933424305;
    coefBTrue[nCoef -  4] =   9.833574586e-05;
    coefBTrue[nCoef -  3] =   6.555706932e-05;
    coefBTrue[nCoef -  2] =   1.638929098e-05;
    coefBTrue[nCoef -  1] =   4.097328656e-06;

    // calculate coefficients using reference and optimised implementation:
    std::vector<real> coefB = gdd.setupCoefB(sample);

    // check right number of coefficients:
    ASSERT_EQ(coefBTrue.size(), coefB.size());

    // assert that correct B coefficients have been computed:
    for(size_t i = 0; i < coefB.size(); i++)
    {
//        std::cout<<"B_true = "<<coefBTrue[i]<<"  "
//                 <<"coefB = "<<coefB[i]<<"  "
//                 <<std::endl;

//        ASSERT_NEAR(coefBTrue[i], coefB[i], tolerance);
    }
}


/*!
 *
 */
TEST_F(GaussianDensityDerivativeTest, GaussianDensityDerivativeConsistencyTest)
{

    std::vector<real> sample = {0.0, 0.33, 0.5, 0.7, 0.4, 0.5, 0.121, 1.0, 2.0};
//    sample = testData_;
    std::vector<real> eval = sample;



    auto minIt = std::min_element(sample.begin(), sample.end());
    auto maxIt = std::max_element(sample.begin(), sample.end());


    real scale = 1.0/((*maxIt) - (*minIt));
    real shift = -(*minIt);
    std::cout<<"max = "<<(*maxIt)<<std::endl;
    std::cout<<"min = "<<(*minIt)<<std::endl;

    std::cout<<"scale = "<<scale<<std::endl;
    std::cout<<"shift = "<<shift<<std::endl;

    std::for_each(sample.begin(), sample.end(), [scale, shift](real &s){s = (s + shift)*scale;});

    for(auto s : sample)
    {
        std::cout<<s<<std::endl;
    }


    real eps = 0.001;

    
    GaussianDensityDerivative gdd;
    gdd.setDerivOrder(2);
    gdd.setBandWidth(1.0);
    gdd.setErrorBound(eps);





    // estimate derivative via direct loop and via approximate method:
    std::vector<real> derivDirect = gdd.estimateDirect(sample, eval);
    std::vector<real> derivApprox = gdd.estimateApprox(sample, eval);
    std::vector<real> deriv = gdd.Evaluate(eval, sample);
    std::vector<real> derivDir = gdd.EvaluateDirect(sample, eval);


    gdd.choose_parameters(sample, eval);


    // loop over all eval points and assert equality of estimation methods:
    for(int i = 0; i < eval.size(); i++)
    {
        std::cout<<"derivDirect = "<<derivDirect[i]<<"  "
                 <<"derivApprox = "<<derivApprox[i]<<"  "
                 <<"deriv = "<<deriv[i]<<"  "
                 <<"derivDir = "<<derivDir[i]<<"  "
                 <<std::endl;
        ASSERT_NEAR(derivDir[i], deriv[i], eps);
        ASSERT_NEAR(derivDir[i], derivDirect[i], eps);
        ASSERT_NEAR(derivDir[i], derivApprox[i], eps);
    }
}

