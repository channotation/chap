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


#ifndef HISTOGRAM_DENSITY_ESTIMATOR_HPP
#define HISTOGRAM_DENSITY_ESTIMATOR_HPP

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

#include "statistics/abstract_density_estimator.hpp"


/*!
 * \brief One-dimensional density estimate using a histogram.
 *
 * This class provides a method estimate() for estimating a probability 
 * density of from a sample of scalar values. It uses a histogram with
 * uniformly spaced bins, where the bin width can be set using setBinWidth().
 *
 * The resulting density is then interpolated linearly using 
 * LinearSplineInterp1D and the result is returned as a SplineCurve1D object.
 * Linear rather than higher order interpolation is used to avoid artifacts
 * where the interpolation method introduces locally negative densities.
 */
class HistogramDensityEstimator : public AbstractDensityEstimator
{
    friend class HistogramDensityEstimatorTest;
    FRIEND_TEST(HistogramDensityEstimatorTest,
                HistogramDensityEstimatorBreaksTest);
    FRIEND_TEST(HistogramDensityEstimatorTest,
                HistogramDensityEstimatorEmptyDatasetTest);
    FRIEND_TEST(HistogramDensityEstimatorTest,
                HistogramDensityEstimatorDensityTest);
    FRIEND_TEST(HistogramDensityEstimatorTest,
                HistogramDensityEstimatorEstimateTest);

    public:
       
        // constructor and destructor:
        HistogramDensityEstimator();

        // density estimation interface:
        virtual SplineCurve1D estimate(
                std::vector<real> &samples);

        // setter methods for parameters:
        virtual void setParameters(
                const DensityEstimationParameters &params);
        void setBinWidth(real binWidth);

    private:

        // internal parameters:
        real binWidth_;

        // auxiliary functions:
        std::vector<real> createBreaks(
                const std::vector<real> &samples);
        std::vector<real> createMidpoints(
                const std::vector<real> &breaks);
        std::vector<real> calculateDensity(
                const std::vector<real> &samples,
                const std::vector<real> &breaks);
    
};

#endif

