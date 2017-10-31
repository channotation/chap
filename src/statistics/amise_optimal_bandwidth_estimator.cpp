#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/tools/roots.hpp>

#include "statistics/amise_optimal_bandwidth_estimator.hpp"
#include "statistics/summary_statistics.hpp"


/*!
 * Estimates the AMISE-optimal bandwidth for kernel density estimation on a 
 * given sample. Requires there to be at least two distinct sample points.
 */
real
AmiseOptimalBandwidthEstimator::estimate(
        const std::vector<real> &samples)
{
    // sanity checks:
    if( samples.size() < 2 )
    {
        throw std::logic_error("Can not use AMISE optimal bandwidth "
                               "estimation with fewer than two samples.");
    }

    // estimate standard deviation of the sample:
    SummaryStatistics sumStats;
    for(auto s : samples)
    {
        sumStats.update(s);
    }
    real sigma = sumStats.sd();

    // estimate functionals six and eight by assuming Gaussianity:
    real phi6 = functionalPhi6(sigma);
    real phi8 = functionalPhi8(sigma);

    // calculate prototype bandwidth optimal wrt asymptotic MSE:
    real g1 = std::pow(-6.0 / ( SQRT2PI_*phi6*samples.size() ), 1.0/7.0);
    real g2 = std::pow(30.0 / ( SQRT2PI_*phi8*samples.size() ), 1.0/9.0);

    // estimate functionals via KDE with bandwidth from optimal MSE:
    real phi4 = functionalPhi(samples, g1, 4);
    phi6 = functionalPhi(samples, g2, 6);

    // calculate constant prefactor in gamma expression:
    gammaFactor_ = gammaFactor(phi4, phi6);

    // parameters for boost root finder:
    real guess = 1.0;
    real factor = 2.0;
    boost::uintmax_t it = 20;
    boost::math::tools::eps_tolerance<real> tol(std::numeric_limits<real>::digits);

    // objective function for root finding:
    std::function<real(real)> objectiveFunction = std::bind(
            &AmiseOptimalBandwidthEstimator::optimalBandwidthEquation,
            this,
            std::placeholders::_1,
            samples);

    // find root:
    std::pair<real, real> root = bracket_and_solve_root(
        objectiveFunction,
        guess, 
        factor, 
        false, 
        tol, 
        it); 

    // return AMISE-optimal bandwidth:
    return root.first;
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
AmiseOptimalBandwidthEstimator::functionalPhi6(real sigma)
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
AmiseOptimalBandwidthEstimator::functionalPhi8(real sigma)
{
    return 105.0 / ( 32.0 * std::pow(sigma, 9) * SQRTPI_ );
}


/*!
 *
 */
real
AmiseOptimalBandwidthEstimator::functionalPhi(
        const std::vector<real> &samples,
        real bw,
        const int deriv)
{
    // initialise sum as zero:
    real phi = 0.0;

    // loop over samples:
    for(auto si : samples)
    {
        for(auto sj : samples)
        {
            // evaluate ermite polynomial:
            real h = boost::math::hermite(deriv, (si - sj)/bw);

            // evaluate kernel itself:
            real k = std::exp(-(si-sj)*(si-sj)/(2.0*bw*bw));

            // add to density derivative functional:
            phi += k*h;
        }
    }

    // multiply constant pre-factor:
    int n = samples.size();
    phi /= n*(n - 1) * SQRT2PI_ * std::pow(bw, deriv + 1);

    // return the density functional:
    return phi;
}


/*!
 * Computes the constant prefactor in gamma(), which needs to be assigned to
 * gammaFactor_ prior to calling gamma().
 */
real
AmiseOptimalBandwidthEstimator::gammaFactor(
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
 * of AmiseOptimalBandwidthEstimator (this precompution must be carried out
 * manually!).
 */
real
AmiseOptimalBandwidthEstimator::gamma(real bw)
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
AmiseOptimalBandwidthEstimator::optimalBandwidthEquation(
        const real bw,
        const std::vector<real> &samples)
{
    // estimate density derivative functional:
    real phi4 = functionalPhi(samples, gamma(bw), 4);

    // return value of implcit equation:
    return bw - std::pow(1.0/( 2.0*SQRTPI_*phi4*samples.size() ), 1.0/5.0);
}

