#include <algorithm>
#include <cmath>

#include <boost/math/special_functions/hermite.hpp>

#include "statistics/gaussian_density_derivative.hpp"

#include <iostream> // TODO


/*!
 * Estimates the density derivative at a set of evaluation points by direct
 * evaluation of the hermite polynomial.
 */
std::vector<real>
GaussianDensityDerivative::estimateDirect(
        const std::vector<real> &sample,
        const std::vector<real> &eval)
{
    // allocate output vector of all zeros:
    std::vector<real> deriv;
    deriv.reserve(eval.size());

    // loop over target points:
    for(auto e : eval)
    {
        deriv.push_back(estimDirectAt(sample, e));
    }

    // scale with correct prefactor:
    real scale = std::pow(-1, r_)/(SQRT2PI_*sample.size()*std::pow(bw_, r_+1));
    std::for_each(deriv.begin(), deriv.end(), [scale](real &d){d *= scale;}); 

    // return vector of density derivative at each evaluation point:
    return deriv;
}


/*!
 * Estimates the density derivative at a single point by direct evaluation.
 */
real
GaussianDensityDerivative::estimDirectAt(
        const std::vector<real> &sample,
        real eval)
{
    real d = 0.0;
    for(auto s : sample)
    {
        real diff = (eval - s)/bw_;
        d += std::exp(-0.5*diff*diff)*boost::math::hermite(r_, diff);
    }
    return d;
}


/*
 *
 */
std::vector<real>
GaussianDensityDerivative::estimateApprox(
        std::vector<real> &sample,
        std::vector<real> &eval)
{
    // shift data and evaluation points:
    real shift = std::min(*std::min_element(sample.begin(), sample.end()),
                          *std::min_element(sample.begin(), sample.end()));
    std::for_each(sample.begin(), sample.end(), [shift](real &s){s -= shift;});
    std::for_each(eval.begin(), eval.end(), [shift](real &e){e -= shift;});

    // scale data, evaluation points and bandwidth:
    real scale = 1.0/std::max(*std::max_element(sample.begin(), sample.end()),
                              *std::max_element(sample.begin(), sample.end()));
    std::for_each(sample.begin(), sample.end(), [scale](real &s){s *= scale;});
    std::for_each(eval.begin(), eval.end(), [scale](real &e){e *= scale;});
    bw_ *= scale;

    // calculate space partitioning (this is data dependent, bc bw_ is scaled):
    centres_ = setupClusterCentres();
    idx_ = setupClusterIndices(sample);

    // compute data dependent coefficients:
    q_ = setupCoefQ(sample.size());
    epsPrime_ = setupScaledTolerance(sample.size());
    rc_ = setupCutoffRadius();
    trunc_ = setupTruncationNumber();
    coefB_ = setupCoefB();

    // allocate output vector of all zeros:
    std::vector<real> deriv;
    deriv.reserve(eval.size());

    // loop over target points:
    for(auto e : eval)
    {
        deriv.push_back(estimApproxAt(sample, e));
    }

    // scale derivative with correct prefactor:
    real fac = q_;
    std::for_each(deriv.begin(), deriv.end(), [fac](real &d){d *= fac;}); 

    // scale back data and evaluation point:
    scale = 1.0/scale;
    std::for_each(sample.begin(), sample.end(), [scale](real &s){ s *= scale; });
    std::for_each(eval.begin(), eval.end(), [scale](real &e){ e *= scale; });
    std::for_each(sample.begin(), sample.end(), [shift](real &s){s += shift;});
    std::for_each(eval.begin(), eval.end(), [shift](real &e){e += shift; });
    bw_ *= scale;

    // return vector of density derivative at each evaluation point:
    return deriv;
}


/*!
 *
 */
real
GaussianDensityDerivative::estimApproxAt(
        const std::vector<real> &sample,
        real eval)
{
    unsigned int sMax = std::floor(static_cast<real>(r_)/2.0);

    //
    real sum = 0.0;
    for(int l = 0; l < 1; l++)
    {
        // sum up to truncation number terms:
        for(unsigned int k = 0; k < trunc_ - 1; k++)
        {
            unsigned int idxA = 0;

            //
            for(unsigned int s = 0; s <= sMax; s++)
            {
                for(unsigned int t = 0; t <= r_ - 2*s; t++)
                {
                    real centre;
                    real d = (eval - centre)/bw_;
                    sum += coefA_.at(idxA) 
                         * coefB_[0] 
                         * std::exp(-0.5*d*d) * std::pow(d, k + r_ - 2*s - t);

                    // increment indeces:
                    idxA++;
                }
            }
        }
    }

    return sum;
}


/*!
 * Sets bandwidth \f$ h > 0 \f$.
 */
void
GaussianDensityDerivative::setBandWidth(real bw)
{
    bw_ = bw;
}


/*!
 * Sets derivative order \f$ r>0 \f$. Also automatically updated the factorial 
 * of \f$ r \f$ and all coefficients that do not also depend on the data.
 */
void
GaussianDensityDerivative::setDerivOrder(unsigned int r)
{
    r_ = r;
    rFac_ = factorial(r);
    coefA_ = setupCoefA();
}


/*!
 * Sets error bound for approximate method \f$ \epsilon > 0 \f$.
 */
void
GaussianDensityDerivative::setErrorBound(real eps)
{
    eps_ = eps;
}


/*!
 * Sets up a vector of equidistent cluster centres covering the unit interval.
 * The cluster spacing is alf the bandwidth.
 */
std::vector<real>
GaussianDensityDerivative::setupClusterCentres()
{
    // minimum number of intervals to cover unit interval:
    ri_ = bw_/2.0;
    numIntervals_ = static_cast<unsigned int>(std::ceil(1.0/ri_));
    ri_ = 1.0/numIntervals_;

    std::vector<real> centres;
    for(int i = 0; i < numIntervals_; i++)
    {
        centres.push_back(i*ri_ + ri_/2.0);
    }

    return centres;
}


/*!
 * Returns a vector of indices into the vector of cluster centres, where each
 * sample is associated with the closest cluster centre.
 */
std::vector<unsigned int>
GaussianDensityDerivative::setupClusterIndices(
        const std::vector<real> &sample)
{
    std::vector<unsigned int> idx;
    idx.reserve(sample.size());
    for(auto s : sample)
    {
        idx.push_back(std::min(static_cast<unsigned int>(std::floor(s/ri_)),
                               numIntervals_ - 1));
    }
    return idx;
}


/*
 *
 */
std::vector<real>
GaussianDensityDerivative::setupCoefA()
{
    // (inclusive) maximum of indices to coefficient matrix:
    unsigned int sMax = std::floor(static_cast<real>(r_)/2.0);
    unsigned int tMax = r_;
    unsigned int sNum = sMax + 1;
    unsigned int tNum = tMax + 1;

    // precompute constant factor in s index:
    std::vector<real> sConstFac;
    sConstFac.reserve(sNum);
    for(int s = 0; s <= sMax; s++)
    {
        sConstFac.push_back(std::pow(2, s) * factorial(s));
    }

    // precompute constant factor in t index:
    std::vector<real> tConstFac;
    tConstFac.reserve(tNum);
    for(int t = 0; t <= tMax; t++)
    {
        tConstFac.push_back(factorial(t));
    }

    // compute coefficients:
    std::vector<real> coefA;
    for(int s = 0; s <= sMax; s++)
    {
        for(int t = 0; t <= r_ - 2*s; t++)
        {
            coefA.push_back( 
                    std::pow(-1, s + t) * rFac_
                    /(sConstFac[s]*tConstFac[t]*factorial(r_ - 2*s - t)));
        }
    }

    return coefA;
}


/*
 *
 */
std::vector<real>
GaussianDensityDerivative::setupCoefB(
        )
{
        
}


/*
 *
 */
real
GaussianDensityDerivative::setupCoefQ(unsigned int n)
{
    return std::pow(-1, r_) / (std::sqrt(2.0*M_PI)*n*std::pow(bw_, r_ + 1));
}


/*
 *
 */
real GaussianDensityDerivative::setupCutoffRadius()
{
    return bw_/2.0 + 2*bw_*std::sqrt(std::log(std::sqrt(rFac_)/epsPrime_));
}


/*
 *
 */
real
GaussianDensityDerivative::setupScaledTolerance(unsigned int n)
{
    return eps_/(std::abs(q_)*n);
}


/*
 *
 */
unsigned int
GaussianDensityDerivative::setupTruncationNumber()
{
    // hardcoded limit for truncation number:
    unsigned int truncMax = 100;

    // factors constant in the loop:
    real bwSq = bw_*bw_;
    real riSq = ri_*ri_;

    // find lowest truncation number that guarantees error bound:
    for(int p = 0; p <= truncMax; p++)
    {
        // calculate error for given truncation number?
        real b = std::min<real>(rc_, ri_ + 0.5*std::sqrt(riSq + 8.0*p*bwSq));
        real d = ri_ - b;
        real err = rFac_/factorial(p) 
                 * std::pow(ri_*b/bwSq, p) 
                 * std::exp(-0.25*d*d/(bwSq));

        // has error bound been reached?
        if( err < epsPrime_ )
        {
            return p;
        }
    }

    // throw exception if error can not be reduced to within desired bound:
    throw std::runtime_error("Could not");
    return -1;
}


/*!
 * Returns the factorial of the given integer.
 */
unsigned int
GaussianDensityDerivative::factorial(unsigned int n)
{
    unsigned int f = 1;
    for(int i = 1; i <= n; i++)
    {
        f *= i;
    }
    return f;
}













