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


#include <limits>

#include "geometry/linear_spline_interp_1D.hpp"
#include "statistics/weighted_kernel_density_estimator.hpp"


/*!
 * Public interface for computing a continuous smooth spline curve that 
 * approximates given discrete function values.
 */
SplineCurve1D
WeightedKernelDensityEstimator::estimate(
        std::vector<real> &samples,
        std::vector<real> &weights)
{    
    // sanity check:
    if( samples.size() != weights.size() )
    {
        throw std::logic_error("Need same number of samples and weights!");
    }

    // create evaluation points:
    std::vector<real> evalPoints = createEvaluationPoints(
            samples);

    // determine weighted density:
    std::vector<real> weightedDensity = calculateWeightedDensity(
            samples,
            weights,
            evalPoints);

    // return weighted density as a spline curve:
    LinearSplineInterp1D Interp;
    return Interp(evalPoints, weightedDensity);
}


/*!
 * Internal evaluation function that computes the Nadaraya-Watson estimate
 * of the smoothing function to the given data points.
 */
std::vector<real>
WeightedKernelDensityEstimator::calculateWeightedDensity(
        std::vector<real> &samples,
        std::vector<real> &weights,
        std::vector<real> &evalPoints)
{
    // set up the density kernel:
    KernelFunctionPointer kernel = KernelFunctionFactory::create(
            kernelFunction_);

    // allocate the density vector:
    std::vector<real> density(evalPoints.size(), 0.0);
    std::vector<real> weightedDensity(evalPoints.size(), 0.0);

    // loop over evaluation points:
    for(size_t i = 0; i < evalPoints.size(); i++)
    {
        // density is sum over kernel distances:
        for(size_t j = 0; j < samples.size(); j++)
        {
            // evaluate kernel function:
            real kern =  kernel -> operator()( 
                    (evalPoints[i] - samples[j])/bandWidth_ );

            // for weighted and unweighted sums:
            density[i] += kern;
            weightedDensity[i] += kern*weights[j];
        }

        // fend of NaNs occuring if density is too close to zero:
        if( density[i] >= std::numeric_limits<real>::epsilon() )
        {
            // Nadaraya-Watson estimate of local function value:
            weightedDensity[i] /= density[i];
        }
    }

    // return density:
    return(weightedDensity);
}

