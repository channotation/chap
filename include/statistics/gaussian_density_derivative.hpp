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


#ifndef GAUSSIAN_DENSITY_DERIVATIVE_HPP
#define GAUSSIAN_DENSITY_DERUVATIVE_HPP

#include <vector>

#include <gtest/gtest.h>

#include <gromacs/utility/real.h>


/*!
 * \brief Calculates derivative of Gaussian kernel density for a given sample.
 *
 * This class provides the two methods estimateDirect() and estimateApprox()
 * which both return the value of the derivative of a nonparametric probability
 * density estimated via a Gaussian kernel. The difference between the two 
 * methods is that the direct approach returns an exact result, but is of 
 * complexity \f$ \mathcal{O}(NM) \f$, whereas the other method returns an
 * approximate derivative with complexity \f$ \mathcal{O}(N + M) \f$, where
 * \f$ N \f$ and \f$ M \f$ are the number of sample and evaluation points 
 * respectively.
 *
 * Both methods require the kernel bandwidth and the derivative order to be 
 * provided with setBandWidth() and setDerivOrder(). In addition, the 
 * approximate method requires and error bound set with setErrorBound(), which
 * will ensure the maximum difference between the approximate and direct 
 * methods, but also increases computational effort for the approximate 
 * evaluation. The approximate method also requires the data, evaluation points
 * and bandwidth to be scaled and shifted such that they lie in the unit 
 * interval and the convenience functions getShiftAndScaleParams(), 
 * shiftAndScale() and shiftAndScaleInverse() are provided as well.
 *
 * The theory underlying the approximate method is explained in the papers
 * "Fast Computation of Kernel Estimators" by Raykar et. al. and "Very Fast
 * Optimal Bandwidth Selection for Univariate Kernel Density Estimation" by
 * Raykar and Duraiswami.
 *
 * \note This class uses double precision internally to avoid erroneous results
 * due to floating point overflow/underflow and the accumulation of rounding
 * errors.
 */
class GaussianDensityDerivative
{
    friend class GaussianDensityDerivativeTest;
    FRIEND_TEST(GaussianDensityDerivativeTest, 
                GaussianDensityDerivativeShiftScaleTest);
    FRIEND_TEST(GaussianDensityDerivativeTest, 
                GaussianDensityDerivativeSpacePartitioningTest);
    FRIEND_TEST(GaussianDensityDerivativeTest, 
                GaussianDensityDerivativeTruncationTest);
    FRIEND_TEST(GaussianDensityDerivativeTest, 
                GaussianDensityDerivativeCoefATest);
    FRIEND_TEST(GaussianDensityDerivativeTest, 
                GaussianDensityDerivativeCoefBTest);
    FRIEND_TEST(GaussianDensityDerivativeTest, 
                GaussianDensityDerivativeConsistencyTest);

    public:

        // public interface for evaluation:
        std::vector<real> estimateApprox(
                const std::vector<real> &sample,
                const std::vector<real> &eval);
        std::vector<real> estimateDirect(
                const std::vector<real> &sample,
                const std::vector<real> &eval);

        // convenience functions for data preparation:
        std::pair<real, real> getShiftAndScaleParams(
                const std::vector<real> &sample,
                const std::vector<real> &eval);
        void shiftAndScale(
                std::vector<real> &vec, 
                real shift, 
                real scale);
        void shiftAndScaleInverse(
                std::vector<real> &vec, 
                real shift, 
                real scale);

        // setter functions:
        void setBandWidth(real bw);
        void setDerivOrder(unsigned int r);
        void setErrorBound(real eps);

    private:

        // internal variables:
        unsigned int numIntervals_;
        unsigned int r_;
        unsigned int rFac_;
        unsigned int trunc_;

        real bw_;
        real eps_;
        real epsPrime_;
        real q_;
        real ri_;
        real rc_;

        std::vector<real> centres_;
        std::vector<real> coefA_;
        std::vector<real> coefB_;
        std::vector<unsigned int> idx_;

        // estimation at an individual evaluation point: 
        real estimDirectAt(
                const std::vector<real> &sample,
                real eval);
        real estimApproxAt(
                real eval);

        // space partitioning:
        std::vector<real> setupClusterCentres();
        std::vector<unsigned int> setupClusterIndices(
                const std::vector<real> &sample);

        // calculation of coefficients:
        std::vector<real> setupCoefA();
        std::vector<real> setupCoefB(const std::vector<real> &sample);
        real setupCoefQ(unsigned int n);
        real setupCutoffRadius();
        real setupScaledTolerance(unsigned int n);
        unsigned int setupTruncationNumber();

        // internal utilities:
        double hermite(
                double x, 
                unsigned int r);
        double factorial(
                double n);
};

#endif

