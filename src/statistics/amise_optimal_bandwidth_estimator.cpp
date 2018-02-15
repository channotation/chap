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


#include <numeric>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/tools/roots.hpp>

#include "statistics/amise_optimal_bandwidth_estimator.hpp"
#include "statistics/summary_statistics.hpp"


/*!
 * Estimates the AMISE-optimal bandwidth for kernel density estimation on a 
 * given sample. Requires there to be at least two distinct sample points.
 */
real
AmiseOptimalBandWidthEstimator::estimate(
        const std::vector<real> &sampleIn)
{
    // make copy of sample:
    // (not so memory efficient, but avoids rescaling the data back)
    auto sample = sampleIn;

    // sanity checks:
    if( sample.size() < 2 )
    {
        // one angstrom returned in this case as default:
        return 0.1;
    }

    // shift and scale data:
    auto ss = gdd_.getShiftAndScaleParams(sample, sample);
    gdd_.shiftAndScale(sample, ss.first, ss.second);

    // estimate standard deviation of the sample:
    SummaryStatistics sumStats;
    for(auto s : sample)
    {
        sumStats.update(s);
    }
    real sigma = sumStats.sd();

    // estimate functionals six and eight by assuming Gaussianity:
    real phi6 = functionalPhi6(sigma);
    real phi8 = functionalPhi8(sigma);

    // calculate prototype bandwidth optimal wrt asymptotic MSE:
    real g1 = std::pow(-6.0 / ( SQRT2PI_*phi6*sample.size() ), 1.0/7.0);
    real g2 = std::pow(30.0 / ( SQRT2PI_*phi8*sample.size() ), 1.0/9.0);

    // estimate functionals via KDE with bandwidth from optimal MSE:
    real phi4 = functionalPhiFast(sample, g1, 4);
    phi6 = functionalPhiFast(sample, g2, 6);
    
    // calculate constant prefactor in gamma expression:
    gammaFactor_ = gammaFactor(phi4, phi6);

    // parameters for boost root finder:
    // initial guess is Silverman's rule of thumb:
    real guess = 1.06*sigma/std::pow(sample.size(), 1.0/5.0);
    real factor = 2.0;
    boost::uintmax_t it = 20;
    boost::math::tools::eps_tolerance<real> tol(std::numeric_limits<real>::digits - 4);

    // objective function for root finding:
    std::function<real(real)> objectiveFunction = std::bind(
            &AmiseOptimalBandWidthEstimator::optimalBandwidthEquation,
            this,
            std::placeholders::_1,
            sample);

    // find root:
    std::pair<real, real> root = bracket_and_solve_root(
        objectiveFunction,
        guess, 
        factor, 
        true, 
        tol, 
        it); 

    // return AMISE-optimal bandwidth (scaled back to original interval):
    return root.first / ss.second;
}


/*!
 * Calculates the eighth order density derivative functional:
 *
 * \f[
 *      \phi_6 = \frac{-15}{16 \sqrt{\pi} \sigma^7}
 * \f]
 *
 * This is based on an approximation assuming a Gaussian probability density.
 */
real
AmiseOptimalBandWidthEstimator::functionalPhi6(real sigma)
{
    return -15.0 / ( 16.0 * std::pow(sigma, 7) * SQRTPI_ );
}


/*!
 * Calculates the eighth order density derivative functional:
 *
 * \f[
 *      \phi_8 = \frac{105}{32 \sqrt{\pi} \sigma^9}
 * \f]
 *
 * This is based on an approximation assuming a Gaussian probability density.
 */
real
AmiseOptimalBandWidthEstimator::functionalPhi8(real sigma)
{
    return 105.0 / ( 32.0 * std::pow(sigma, 9) * SQRTPI_ );
}


/*!
 * Evaluation the functional
 *
 * \f[
 *      \Phi_r = \frac{1}{N}\sum_{i=1}^N p^{(r)}(s_i)
 * \f]
 *
 * i.e. the sum over the density derivative evaluated at each sample point. The
 * derivative itself is evaluated using the GaussianDensityDerivative class,
 * which uses an approximate method to evaluate this expression in linear 
 * complexity with the number of sample points.
 */
real
AmiseOptimalBandWidthEstimator::functionalPhiFast(
        const std::vector<real> &sample,
        real bw,
        int deriv)
{
    // set density derivative estimation parameters:
    gdd_.setErrorBound(0.01);
    gdd_.setBandWidth(bw);
    gdd_.setDerivOrder(deriv);

    // obtain derivative estimate at each sample point:
    std::vector<real> d = gdd_.estimateApprox(sample, sample);
    
    // average of derivatives is estimate for phi:
    real phi = std::accumulate(d.begin(), d.end(), 0.0);
    phi /= sample.size();

    // return phi parameter:
    return phi;
}


/*!
 * Computes the constant prefactor in gamma(), which needs to be assigned to
 * gammaFactor_ prior to calling gamma().
 */
real
AmiseOptimalBandWidthEstimator::gammaFactor(
        const real phi4,
        const real phi6)
{
    return std::pow(-6.0*std::sqrt(2.0)*phi4 / phi6 , 1.0/7.0);
}


/*!
 * Returns the bandwidth, \f$ \gamma \f$, used to estimate the density 
 * derivative functional entering the optimalBandwidthEquation():
 *
 * \f[
 *      \gamma = \left[ \frac{-6\sqrt{2}\Phi_4(g_1)}{\Phi_6(g_2)} \right]^{\frac{1}{7}} h^{\frac{5}{7}}
 * \f]
 *
 * Note that as the prefactor on square brackets is independent of \f$ h \f$,
 * it is computed only once using gammaFactor() and stored in a member variable
 * of AmiseOptimalBandWidthEstimator (this precompution must be carried out
 * manually!).
 */
real
AmiseOptimalBandWidthEstimator::gamma(real bw)
{
    return gammaFactor_ * std::pow(bw, 5.0/7.0);
}


/*!
 * Returns the value for the implicit expression for the AMISE-optimal
 * bandwidth:
 *
 * \f[
 *      h - \left[ \frac{1.0}{2\sqrt{\pi}\Phi_4\big(\gamma (h)\big) N} \right]^{\frac{1}{5}}
 * \f]
 *
 * This is solved iteratively to obtain \f$ h \f$, which can then used to 
 * obtain a kernel density estimate via the KernelDensityEstimator class.
 */
real
AmiseOptimalBandWidthEstimator::optimalBandwidthEquation(
        const real bw,
        const std::vector<real> &samples)
{

    // estimate density derivative functional:
    real phi4 = functionalPhiFast(samples, gamma(bw), 4);

    // evaluate optimal bandwidth equation:
    return bw - std::pow(1.0/(2.0*SQRTPI_*phi4*samples.size()) , 1.0/5.0);
}

