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


#ifndef KERNEL_DENSITY_ESTIMATOR_HPP
#define KERNEL_DENSITY_ESTIMATOR_HPP

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/real.h"

#include "geometry/spline_curve_1D.hpp"
#include "statistics/abstract_density_estimator.hpp"
#include "statistics/kernel_function.hpp"


/*!
 * \brief One-dimensional density estimation using a kernel formalism.
 *
 * This class provides a method estimate() for estimating the probability 
 * density from a sample of scalar values. It uses one of the classes derived
 * from AbstractKernelFunction to calculate the density according to
 *
 * \f[
 *      p(x_j) = \frac{1}{h N} \sum_{i=1}^{N} K\left( \frac{x_j - x_i}{h} \right)
 * \f]
 *
 * with a band width \f$ h \f$. The evaluation points \f$ x_j \f$ are spaced
 * uniformly and are no further than a user specified distance apart. 
 *
 * The resulting density is interpolated linearly using LinearSplineInterp1D
 * in order to avoid overshoots resulting in negative densities that may 
 * occur with higher order interpolation.
 */
class KernelDensityEstimator : public AbstractDensityEstimator
{
    friend class KernelDensityEstimatorTest;
    FRIEND_TEST(
            KernelDensityEstimatorTest, 
            KernelDensityEstimatorParameterTest);
    FRIEND_TEST(
            KernelDensityEstimatorTest, 
            KernelDensityEstimatorEvalPointTest);
    FRIEND_TEST(
            KernelDensityEstimatorTest, 
            KernelDensityEstimatorGaussianRawDensityTest);
    FRIEND_TEST(
            KernelDensityEstimatorTest, 
            KernelDensityEstimatorEmptyDatasetTest);
    FRIEND_TEST(
            KernelDensityEstimatorTest, 
            KernelDensityEstimatorGaussianInterpDensityTest);

    public:
        
        // density estimation interface:
        virtual SplineCurve1D estimate(
                std::vector<real> &samples);

        // implementation of parameter setting interface:
        virtual void setParameters(
                const DensityEstimationParameters &params);

    protected:

        // flag indicating whether parameters have been set:
        bool paramtersSet_;

        // internal parameters:
        real bandWidth_;
        real bandWidthScale_;
        real maxEvalPointDist_;
        real evalRangeCutoff_;
        eKernelFunction kernelFunction_;

        // auxiliary functions for parameter setting:
        void setBandWidth(const real bandWidth);
        void setBandWidthScale(const real scale);
        void setMaxEvalPointDist(const real maxEvalPointDist);
        void setEvalRangeCutoff(const real evalRangeCutoff);
        void setKernelFunction(const eKernelFunction kernelFunction);

        // auxiliary functions for density estimation:
        std::vector<real> createEvaluationPoints(
                const std::vector<real> &samples);
        size_t calculateNumEvalPoints(
                const real range);
        std::vector<real> calculateDensity(
                const std::vector<real> &samples,
                const std::vector<real> &evalPoints);
        void endpointDensityToZero(
                std::vector<real> &density,
                std::vector<real> &evalPoints);
};

#endif

