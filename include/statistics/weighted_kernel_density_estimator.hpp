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


#ifndef WEIGHTED_KERNEL_DENSITY_ESTIMATOR_HPP
#define WEIGHTED_KERNEL_DENSITY_ESTIMATOR_HPP

#include <vector>

#include "geometry/spline_curve_1D.hpp"
#include "statistics/kernel_density_estimator.hpp"


/*!
 * \brief Nadaraya-Watson kernel smoother.
 *
 * Performs essentially the same job as KernelDensityEstimator, but instead
 * of merely returning an estimate for the density of a point set that is 
 * normalised such that it integrates to one, it returns a continues function
 * that aims at smoothly interpolating between the given function values 
 * (interpreted as weights at the sample values).
 */
class WeightedKernelDensityEstimator : public KernelDensityEstimator
{
    public:

        // estimation interface:
        SplineCurve1D estimate(
                std::vector<real> &samples,
                std::vector<real> &weights);


    private:

        // evaluation of weighted density:
        std::vector<real> calculateWeightedDensity(
                std::vector<real> &samples,
                std::vector<real> &weights,
                std::vector<real> &evalPoints);
};

#endif

