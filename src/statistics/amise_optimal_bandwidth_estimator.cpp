#include <algorithm>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/tools/roots.hpp>

#include "statistics/amise_optimal_bandwidth_estimator.hpp"
#include "statistics/summary_statistics.hpp"

using namespace boost::math;


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
    real phi6 = functionalPhi6(sigma); // TODO correct!
    real phi8 = functionalPhi8(sigma); // TODO correct!

    // calculate prototype bandwidth optimal wrt asymptotic MSE:
    real g1 = std::pow(-6.0 / ( SQRT2PI_*phi6*samples.size() ), 1.0/7.0); // TODO correct!
    real g2 = std::pow(30.0 / ( SQRT2PI_*phi8*samples.size() ), 1.0/9.0); // TODO correct!

/*
    std::cout<<std::endl<<std::endl;
    std::cout<<"phi6 = "<<phi6<<"  "
             <<"phi8 = "<<phi8<<"  "
             <<"g1 = "<<g1<<"  "
             <<"g2 = "<<g2<<"  "
             <<std::endl; */

    // estimate functionals via KDE with bandwidth from optimal MSE:
    real phi4 = functionalPhi(samples, g1, 4); // TODO correct!
    phi6 = functionalPhi(samples, g2, 6); // TODO correct!

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

    // parameters for boost root finder:
    std::cout<<std::endl;
    real guess = 1.0; // TODO!
    real factor = 2.0;
    boost::uintmax_t it = 20;
    boost::math::tools::eps_tolerance<real> tol(std::numeric_limits<real>::digits - 3);

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
        true, 
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



/*
 *
 */
real
AmiseOptimalBandwidthEstimator::functionalPhiFast(
        std::vector<real> &samples,
        real bw,
        const int deriv)
{
    // TODO: handle this as input?
    real eps = 0.01;

    clock_t t_prep = std::clock();

    // sort samples (TODO: move this to higher level function?):
    std::sort(samples.begin(), samples.end());


    // shift and scale the input data (assumes samples to be sorted!):
    real shift = *std::min_element(samples.begin(), samples.end());
    std::for_each(samples.begin(), samples.end(), [shift](real &s){s-=shift;});
    real scale = 1.0 / *std::max_element(samples.begin(), samples.end());
    std::for_each(samples.begin(), samples.end(), [scale](real &s){s*=scale;});

/*
    for(auto s : samples)
    {
        std::cout<<"s = "<<s<<std::endl;
    }*/

    bw = scale*bw;

/*
    std::cout<<"scale = "<<scale<<"  "
             <<"shift = "<<shift<<"  "
             <<"deriv = "<<deriv<<"  "
             <<std::endl;*/

    // calculate q-factor:
    real q = coefQ(bw, deriv, samples.size());
    real epsPrime = eps / (samples.size() * std::abs(q));

    // calculate interval radius:
    real ri = 0.5*bw;

    // find interval centres and cluster point by interval centre:
    std::vector<real> centres = intervalCentres(bw);
    std::vector<size_t> centreIdx = nearestIntervalCentre(centres, samples);
    
    // calculate factorial of derivative:
    real derivFactorial = factorial<real>(deriv);

    // calculate curoff radius:
    real rc = cutoffRadius(ri, bw, derivFactorial, epsPrime);

//    std::cout<<"pre truncation number"<<std::endl;    

    //  calculate truncation number:
    int p = truncationNumber(bw, ri, rc, epsPrime, deriv);

//    std::cout<<"pre approx density derivative"<<std::endl;
    

    t_prep = std::clock() - t_prep;

    clock_t t_phi = std::clock();

    //
    real phi = 0.0;
    for(auto s : samples)
    {
//        std::cout<<"phi = "<<phi<<std::endl;

        phi += approximateDensityDerivative(
                samples,
                centres,
                centreIdx,
                s,
                bw,
                rc,
                ri,
                epsPrime,
                q,
                deriv,
                p);
    }

    int n = samples.size();
    phi /= (n-1);

    t_phi = std::clock() - t_phi;

    std::cout<<std::endl;
    std::cout<<"t_prep = "<<t_prep<<std::endl;
    std::cout<<"t_phi = "<<t_phi<<std::endl;

    // return phi:
    return phi;
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
    real bwIn = bw;

    // estimate density derivative functional:
    real phi4 = functionalPhi(samples, gamma(bw), 4);

    // return value of implcit equation:
//    real val = bw - std::pow(1.0/( 2.0*SQRTPI_*phi4*samples.size() ), 1.0/5.0);

    real val = bw - std::pow(1/(2*SQRTPI_*phi4*samples.size()) , 1.0/5.0);

    std::cout<<"bw = "<<bw<<"  "
             <<"bwIn = "<<bwIn<<"  "
             <<"gammaFac = "<<gammaFactor_<<"  "
             <<"gamma = "<<gamma(bw)<<"  "
             <<"phi4 = "<<phi4<<"  "
             <<"val = "<<val<<std::endl;

    return val;
}


/*!
 * Returns a vector of interval centres covering the unit interval 
 * \f$ [0.0, 1.0] \f$, where the distance between two interval centres is equal
 * to the given bandwidth. The position of the lowest interval centre is 
 * determined by
 *
 * \f[
 *      c_0 = \frac{M*h - 1}{2}
 * \f]
 *
 * where \f$ M \f$ is the miniumum number of intervals required to cover the 
 * entire unit interval.
 */
std::vector<real>
AmiseOptimalBandwidthEstimator::intervalCentres(const real bw)
{
    // how many intervals do we need to giver unit interval?
    real numInterv = std::ceil(1.0/bw);

    // factor in lowest overhand of covering intervals:
    real centreLo = -0.5*(bw*numInterv - 1.0);

    // build vector of interval centres:
    std::vector<real> intervalCentres;
    intervalCentres.reserve(numInterv);
    for(size_t i = 0; i < numInterv; i++)
    {
        intervalCentres.push_back(centreLo + i*bw);
    }

    // return the interval centres:
    return intervalCentres;
}


/*!
 * Finds the index of the nearest centre for each sample. Requires the samples
 * vector to be sorted! Note that the assignment of sample points to clusters 
 * is precomputed using this function rather than assigning cluster on the fly,
 * because this enables a consistent treatment of data points located halfway
 * between two cluster centres.
 */
std::vector<size_t>
AmiseOptimalBandwidthEstimator::nearestIntervalCentre(
        const std::vector<real> &centres,
        const std::vector<real> &samples)
{
    // allocate memory for nearest centre indices:
    std::vector<size_t> centreIdx;
    centreIdx.reserve(samples.size());

    // loop over samples:
    for(auto &s : samples)
    {
        // find an upper bound:
        auto ub = std::lower_bound(centres.begin(), centres.end(), s);

        // handle edge cases:
        if( ub == centres.begin() )
        {
            centreIdx.push_back(0);
            continue;
        }
        if( ub == centres.end() )
        {
            centreIdx.push_back(centres.size() - 1);
            continue;
        }

        // get lower bound (simply previous centre):
        auto lb = ub--;

        // which is closer?
        if( *ub - s < s - *lb )
        {
            centreIdx.push_back(std::distance(ub, centres.begin()));
        }
        else
        {
            centreIdx.push_back(std::distance(lb, centres.begin()));
        }
    }

    // return vector of nearest centre indices:
    return centreIdx;
}


/*!
 * Returns the truncation number that guarantees a given error bound.
 */
int
AmiseOptimalBandwidthEstimator::truncationNumber(
        const real bw,
        const real ir,
        const real cr,
        const real epsPrime,
        const unsigned int deriv)
{
    // hard coded limit to guarantee termination:
    int maxTruncNum = 100;

    // calculate squared bandwidth:
    real sqbw = bw*bw;

    // increment truncation number until error bound is met:
    // TODO: a bisection based method would be faster!
    int p = 1;
    while(p <= maxTruncNum)
    {
//        std::cout<<"p = "<<p<<std::endl;

        // calculate b smaller then cutoff radius:
        real b = std::min<real>(cr, (ir + std::sqrt(ir + 8.0*p*sqbw))/2.0);

        // compute error bound:
        real tmp = std::exp( -(ir-b)*(ir-b)/(4.0*sqbw) );
        tmp *= std::pow(ir*b/sqbw, p);
        tmp *= std::sqrt(deriv)/boost::math::factorial<real>(p);

        // target error bound reached?
        if( tmp <= epsPrime )
        {
            return p;
        }

        // increment trial truncation number:
        p++;
    }

    // if we haven't returned yet, target error could not be met:
    throw std::logic_error("Could not meet truncation error limit!");

    // return negative truncation number:
    return -1;
}


/*!
 * Calculates the cutoff radius according to:
 *
 * \f[
 *      r_{\text{c}} = r_{i} + 2h\sqrt{ \log{\sqrt{r!}}/\epsilon^\prime }
 * \f]
 */
real
AmiseOptimalBandwidthEstimator::cutoffRadius(
        const real &ri,
        const real &bw,
        const real &derivFactorial,
        const real &epsPrime)
{
    return ri + 2.0*bw*std::sqrt(std::log(std::sqrt(derivFactorial)/epsPrime));
}


/*!
 *
 */
real
AmiseOptimalBandwidthEstimator::approximateDensityDerivative(
        const std::vector<real> &samples,
        const std::vector<real> &centres,
        const std::vector<size_t> &centreIdx,
        real eval,
        real bw,
        real rc,
        real ri,
        real epsPrime,
        real q,
        int deriv,
        int p)
{
    int test = 0;
    int total = 0;
    // loop over all centres:
    real densityDerivative = 0.0;
    for(int l = 0; l < centres.size(); l++)
    {

        // only perform addition for eval point within cutoff radius:
        if( std::abs(eval - centres[l]) > rc )
        {
 /*           std::cout<<"  "
                     <<"eval = "<<eval<<"  "
                     <<"centre = "<<centres[l]<<"  "
                     <<"rc = "<<rc<<"  "
                     <<std::endl;*/
            continue;
        }

        test++;

        // sum over Taylor series up to truncation number:
        for(size_t k = 0; k < p - 1; k++)
        {
            for(size_t s = 0; s < std::floor(deriv/2); s++)
            {
                for(size_t t = 0; t < deriv - 2*s; t++)
                {
                    real tmp = coefA(s, t, deriv);
//                    std::cout<<" tmp_a = "<<tmp<<"  ";
                    tmp *= coefB(k, t, l, samples, centres, centreIdx, bw, q);
//                    std::cout<<" tmp_b = "<<tmp<<"  ";
                    tmp *= coefC(k, s, t, deriv, eval, centres[l], bw);
//                    std::cout<<" tmp_c = "<<tmp<<"  "<<std::endl;
                    densityDerivative += tmp;

                    total++;
                }
            }
        }
//        std::cout<<"l = "<<l<<"  "
//                 <<"densityDerivative = "<<densityDerivative<<std::endl;
    }


    std::cout<<"p = "<<p<<"  "
             <<"centres.size = "<<centres.size()<<"  "
             <<"centres.used = "<<test<<"  "
             <<"deriv = "<<deriv<<"  "
             <<"totalComplexity = "<<total<<"  "
             <<"N = "<<samples.size()<<"  "
             <<std::endl;

    // return density derivative:
    return densityDerivative;
}


/*!
 * Returns the following coefficient:
 *
 * \f[
 *      a_{st} = \frac{-1^{s+t} r!}{2^s s! t! (r - 2s -t)!}
 * \f]
 */
real
AmiseOptimalBandwidthEstimator::coefA(
        int s,
        int t,
        int deriv)
{
    real tmp =  std::pow(-1, s + t) * factorial<real>(deriv) / (std::pow(2, s) 
           * factorial<real>(s) * factorial<real>(t) * factorial<real>(deriv - 2*s - t));
/*    std::cout<<"coefA = "<<tmp<<"  "
             <<"s = "<<s<<"  "
             <<"t = "<<t<<"  "
             <<"deriv = "<<deriv<<"  "
             <<"s! = "<<std::tgamma(s)<<"  "
             <<"t! = "<<factorial<real>(t)<<"  "
             <<std::endl;*/
    return tmp;
}


/*!
 * Returns the following coefficient:
 *
 * \f[
 *      B_{kt}^l = \frac{q}{k!} \sum_{x_i\in S_l} \exp\left( -\frac{(x_i - c_l)}{2h^2} \right) \left( \frac{x_i - c_l}{h} \right)^{k+t}
 * \f]
 */
real
AmiseOptimalBandwidthEstimator::coefB(
        int k,
        int t,
        int l,
        const std::vector<real> &samples,
        const std::vector<real> &centres,
        const std::vector<size_t> &centreIdx,
        real bw,
        real coefQ)
{
    // loop over all samples associated with this centre:
    real coefB = 0.0;
    for(size_t i = 0; i < samples.size(); i++)
    {
        // is sample in l-th interval?
        if( centreIdx[i] == l )
        {
            real tmp = (samples[i] - centres[l]) / bw;
            coefB += std::exp( -0.5*tmp*tmp ) * std::pow(tmp, k + t);
        }
    }

    // scale and retrun coefficient:
    coefB *= coefQ/factorial<real>(k);             
    return coefB;
}


/*!
 * Returns the following coefficient:
 *
 * \f[
 *      c_l = \exp\left( -\frac{(x_j - c_l)^2}{2h^2} \right) \left( \frac{x_j - c_l}{h} \right)^{k+r-2s-t}
 * \f]
 */
real
AmiseOptimalBandwidthEstimator::coefC(
        int k,
        int s,
        int t,
        int deriv,
        real eval,
        real centre,
        real bw)
{
    real tmp = (eval - centre)/bw;
    return std::exp(-0.5*tmp*tmp) * std::pow(tmp, k + deriv - 2*s -t);
}


/*!
 * Returns the following coefficient:
 *
 * \f[
 *      q = \frac{(-1)^r}{\sqrt{2\pi}Nh^{r + 1}}
 * \f]
 */
real
AmiseOptimalBandwidthEstimator::coefQ(
        real bw,
        int deriv,
        int numSamples)
{
    return std::pow(-1.0, deriv)/(SQRT2PI_*numSamples*std::pow(bw, deriv + 1));
}



