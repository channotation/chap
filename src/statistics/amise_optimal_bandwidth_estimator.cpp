#include <algorithm>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/tools/roots.hpp>

#include "statistics/amise_optimal_bandwidth_estimator.hpp"
#include "statistics/summary_statistics.hpp"

using namespace boost::math;


#include <chrono> // TODO
using namespace std::chrono;


/*!
 * Estimates the AMISE-optimal bandwidth for kernel density estimation on a 
 * given sample. Requires there to be at least two distinct sample points.
 */
real
AmiseOptimalBandwidthEstimator::estimate(
        const std::vector<real> &sampleIn)
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    // make copy of sample:
    auto sample = sampleIn;

    // sanity checks:
    if( sample.size() < 2 )
    {
        throw std::logic_error("Can not use AMISE optimal bandwidth "
                               "estimation with fewer than two samples.");
    }

    high_resolution_clock::time_point t4 = high_resolution_clock::now();

    // shift and scale data:
    auto ss = gdd_.getShiftAndScaleParams(sample, sample);
    gdd_.shiftAndScale(sample, ss.first, ss.second);

    high_resolution_clock::time_point t5 = high_resolution_clock::now();

    // estimate standard deviation of the sample:
    SummaryStatistics sumStats;
    for(auto s : sample)
    {
        sumStats.update(s);
    }
    real sigma = sumStats.sd();

    high_resolution_clock::time_point t6 = high_resolution_clock::now();

    // estimate functionals six and eight by assuming Gaussianity:
    real phi6 = functionalPhi6(sigma); // TODO correct!
    real phi8 = functionalPhi8(sigma); // TODO correct!

    // calculate prototype bandwidth optimal wrt asymptotic MSE:
    real g1 = std::pow(-6.0 / ( SQRT2PI_*phi6*sample.size() ), 1.0/7.0); // TODO correct!
    real g2 = std::pow(30.0 / ( SQRT2PI_*phi8*sample.size() ), 1.0/9.0); // TODO correct!

/*
    std::cout<<std::endl<<std::endl;
    std::cout<<"phi6 = "<<phi6<<"  "
             <<"phi8 = "<<phi8<<"  "
             <<"g1 = "<<g1<<"  "
             <<"g2 = "<<g2<<"  "
             <<std::endl; */
    high_resolution_clock::time_point t7 = high_resolution_clock::now();

    // estimate functionals via KDE with bandwidth from optimal MSE:
    real phi4 = functionalPhiFast(sample, g1, 4); // TODO correct!
    phi6 = functionalPhiFast(sample, g2, 6); // TODO correct!
    
    high_resolution_clock::time_point t8 = high_resolution_clock::now();

/*
    std::cout<<std::endl<<std::endl;
    std::cout<<"phi6 = "<<phi6<<"  "
             <<"phi8 = "<<phi8<<"  "
             <<"g1 = "<<g1<<"  "
             <<"g2 = "<<g2<<"  "
             <<"phi4 = "<<phi4<<"  "
             <<"phi6 = "<<phi6<<"  "
             <<std::endl;*/

    // calculate constant prefactor in gamma expression:
    gammaFactor_ = gammaFactor(phi4, phi6); // TODO correct!

    high_resolution_clock::time_point t9 = high_resolution_clock::now();
    

    // parameters for boost root finder:
    // initial guess is Silverman's rule of thumb:
    real guess = 1.06*sigma/std::pow(sample.size(), 1.0/5.0);
    real factor = 2.0;
    boost::uintmax_t it = 20;
    boost::math::tools::eps_tolerance<real> tol(2); // two significant sigits

    // objective function for root finding:
    std::function<real(real)> objectiveFunction = std::bind(
            &AmiseOptimalBandwidthEstimator::optimalBandwidthEquation,
            this,
            std::placeholders::_1,
            sample);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    // find root:
    std::pair<real, real> root = bracket_and_solve_root(
        objectiveFunction,
        guess, 
        factor, 
        true, 
        tol, 
        it); 

    high_resolution_clock::time_point t3 = high_resolution_clock::now();

    
    auto prep = duration_cast<microseconds>( t2 - t1 ).count();
    auto opti = duration_cast<microseconds>( t3 - t2 ).count();
    auto tscale = duration_cast<microseconds>( t4 - t1 ).count();
    auto tsumStats = duration_cast<microseconds>( t5 - t4 ).count();
    auto tphi46 = duration_cast<microseconds>( t8 - t7 ).count();
    auto tgammafac = duration_cast<microseconds>( t9 - t8 ).count();
    auto trootprep = duration_cast<microseconds>( t2 - t9 ).count();

    std::cout<<"prep = "<<prep<<"  "
             <<"opti = "<<opti<<"  "
             <<"prop/opti = "<<(double)prep/(double)opti*100<<"  "
             <<"scale = "<<tscale<<"  "
             <<"sumStats = "<<tsumStats<<"  "
             <<"phi46 = "<<tphi46<<"  "
             <<"gamm = "<<tgammafac<<"  "
             <<"rootprep = "<<trootprep<<"  "
             <<std::endl;


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



/*
 *
 */
real
AmiseOptimalBandwidthEstimator::functionalPhiFast(
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
 *
 */
/*
real
AmiseOptimalBandwidthEstimator::functionalPhi(
        const std::vector<real> &samples,
        real bw,
        const int deriv)
{
    // initialise sum as zero:
    double phi = 0.0;

    // loop over samples:
    for(auto si : samples)
    {
        for(auto sj : samples)
        {
            // evaluate hermite polynomial:
            real h = boost::math::hermite(deriv, (si - sj)/bw/std::sqrt(2.0)) 
                   * std::pow(2, -deriv/2);

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
}*/


/*!
 * Computes the constant prefactor in gamma(), which needs to be assigned to
 * gammaFactor_ prior to calling gamma().
 */
real
AmiseOptimalBandwidthEstimator::gammaFactor(
        const real phi4,
        const real phi6)
{
//    std::cout<<"phi4 = "<<phi4<<"  "
//             <<"phi6 = "<<phi6<<"  "
//             <<std::endl;
             // FIXME ERROE here!
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
    real phi4 = functionalPhiFast(samples, gamma(bw), 4);

    // return value of implcit equation:
//    real val = bw - std::pow(1.0/( 2.0*SQRTPI_*phi4*samples.size() ), 1.0/5.0);

    real val = bw - std::pow(1.0/(2.0*SQRTPI_*phi4*samples.size()) , 1.0/5.0);


    std::cout<<"bw = "<<bw<<"  "
             <<"gammaFac = "<<gammaFactor_<<"  "
             <<"gamma = "<<gamma(bw)<<"  "
             <<"phi4 = "<<phi4<<"  "
             <<"val = "<<val<<std::endl;

    return val;
}


