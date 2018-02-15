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


#include <algorithm>
#include <cmath>

#include "statistics/gaussian_density_derivative.hpp"


/*!
 * Estimates the density derivative at a set of evaluation points by direct
 * evaluation of the hermite polynomial. This has an overall complexity 
 * propertional to the product of number of samples and number of evaluation
 * points and is hence of quadratic complexity if evalution points and sample
 * points are identical.
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
    real scale = setupCoefQ(sample.size());
    std::for_each(deriv.begin(), deriv.end(), [scale](real &d){d *= scale;}); 

    // return vector of density derivative at each evaluation point:
    return deriv;
}


/*!
 * Estimates the density derivative at a single point by direct evaluation
 * of
 *
 * \f[
 *      p^{(r)}(e) \frac{(-1)^r}{\sqrt{2\pi}nh^{r+1}} \sum_{i=1}^n H_r\left( \frac{e - s_i}{h} \right) \exp\left( -\frac{(e - s_i)^2}{2h^2} \right)
 * \f]
 *
 * where \f$ H_r(x) \f$ is the probabilist's hermite polynomial.
 */
real
GaussianDensityDerivative::estimDirectAt(
        const std::vector<real> &sample,
        real eval)
{
    double d = 0.0;
    for(auto s : sample)
    {
        double diff = (eval - s)/bw_;
        d += std::exp(-0.5*diff*diff)*hermite(diff, r_);
    }
    return d;
}


/*!
 * Evaluates the derivative of a Gaussian kernel density using an approximate 
 * expression that is of linear complexity in the order of samples and
 * evaluation points. This function calculates coefficients and then passes
 * the evaluation of the derivative at one given evaluation point to 
 * estimateApproxAt(). Note that this function assumes the data to lie in the 
 * interval \f$ [0,1] \f$.
 */
std::vector<real>
GaussianDensityDerivative::estimateApprox(
        const std::vector<real> &sample,
        const std::vector<real> &eval)
{
    // calculate space partitioning (this is data dependent, bc bw_ is scaled):
    centres_ = setupClusterCentres();
    idx_ = setupClusterIndices(sample);

    // compute data dependent coefficients:
    q_ = setupCoefQ(sample.size());
    epsPrime_ = setupScaledTolerance(sample.size());
    rc_ = setupCutoffRadius();
    trunc_ = setupTruncationNumber();
    coefB_ = setupCoefB(sample);

    // loop over target points:
    std::vector<real> deriv;
    
    deriv.reserve(eval.size());
    for(auto e : eval)
    {
        deriv.push_back(estimApproxAt(e));
    }
    
    // return vector of density derivative at each evaluation point:
    return deriv;
}


/*!
 * Evaluates the \f$ r \f$-th derivative of the Gaussian density at the given
 * evaluation point using the approximate expression
 *
 * \f[
 *      p^{(r)}_\epsilon(e) \sum_{ l : \left| e - c_l \right| \leq r_\text{c} } \sum_{k=0}^{p-1} \sum_{s=0}^{ \lfloor r/2 \rfloor } \sum_{t=0}^{r-2s} a_{st} B_{kt}^l \exp\left( -\frac{(e - c_l)^2}{2h^2} \right) \times \left( \frac{e - c_l}{h} \right)^{k + r - 2s - t}
 * \f]
 *
 * where the coefficients \f$ a_{st} \f$ and \f$ B_{kt}^l \f$ can be 
 * precomputed for repeated use of this function.
 */
real
GaussianDensityDerivative::estimApproxAt(
        real eval)
{ 
    // upper bound for coefficient loop:
    unsigned int sMax = floor(static_cast<real>(r_)/2.0);

    // sum up the terms in approximation of derivative:
    double sumPos = 0.0;
    double sumNeg = 0.0;
    double sum = 0.0;           
    for(unsigned int l = 0; l < centres_.size(); l++)
    {
        // distance from cluster centre:
        double dist = eval - centres_[l];

        // skip this iteration if cluster centre beyond cutoff radius:
        if( std::fabs(dist) > rc_ )
        {
            continue;
        }

        // scale distance with bandwidth:
        dist /= bw_;

        // precompute exponential term:
        double expTerm = exp(-0.5*dist*dist);

        // also precompute power term:
        std::vector<double> powTerm(trunc_ + r_);
        powTerm[0] = 1.0;
        for(unsigned int i = 1; i < trunc_ + r_; i++)
        {   
            powTerm[i] = powTerm[i-1] * dist;
        }   

        // loop up to truncation number:
        for(unsigned int k = 0; k <= trunc_ - 1; k++)
        {   
            // A-coefficients will be accessed in order of creation:
            unsigned int idxA = 0;

            // loop over coefficients:
            for(unsigned int s = 0; s <= sMax; s++)
            {   
                for(unsigned int t = 0; t <= r_ - 2*s; t++)
                {   
                    // add to sum:
                    sum += coefA_[idxA]
                         * coefB_[l*trunc_*(r_ + 1) + (r_ + 1)*k + t]
                         * expTerm
                         * powTerm[k + r_ - 2*s - t];

                    // increment A-coefficient index:
                    idxA++;
                }   
            }   
        }   
    }   

    // return density derivative at evaluation point:
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
 * Sets up a vector of equidistant cluster centres covering the unit interval.
 * The cluster spacing is half the bandwidth.
 */
std::vector<real>
GaussianDensityDerivative::setupClusterCentres()
{
    // minimum number of intervals to cover unit interval:
    ri_ = bw_/2.0;
    numIntervals_ = static_cast<unsigned int>(std::ceil(1.0/ri_));
    ri_ = 1.0/numIntervals_;

    std::vector<real> centres;
    for(unsigned int i = 0; i < numIntervals_; i++)
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


/*!
 * Calculates the coefficients
 *
 * \f[
 *      a_{st} = \frac{(-1)^{s+t} r!}{2^2 s! t! (r - 2s - t)!}
 * \f]
 *
 * which depend only on the derivative order and can be reused with new data.
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
    for(unsigned int s = 0; s <= sMax; s++)
    {
        sConstFac.push_back(std::pow(2, s) * factorial(s));
    }

    // precompute constant factor in t index:
    std::vector<real> tConstFac;
    tConstFac.reserve(tNum);
    for(unsigned int t = 0; t <= tMax; t++)
    {
        tConstFac.push_back(factorial(t));
    }

    // compute coefficients:
    std::vector<real> coefA;
    for(unsigned int s = 0; s <= sMax; s++)
    {
        for(unsigned int t = 0; t <= r_ - 2*s; t++)
        {
            coefA.push_back( 
                    std::pow(-1, s + t) * rFac_
                    /(sConstFac[s]*tConstFac[t]*factorial(r_ - 2*s - t)));
        }
    }

    return coefA;
}


/*!
 * Calculates the coefficient matrix
 *
 * \f[
 *      B_{kt}^l = \sum_{s_i \in S_l} q \exp\left( -\frac{(x_i - c_l)^2}{2h^2} \right) \times  \left( \frac{s_i - c_l}{h} \right)^{k+1}
 * \f]
 *
 * where \f$ S_l \f$ is the \f$ l \f$-th interval computed with 
 * setupClusterCentres(). Note that these coefficients have to be recomputed 
 * for new data.
 */
std::vector<real>
GaussianDensityDerivative::setupCoefB(
        const std::vector<real> &sample)
{
    // allocate coefficient matrix:
    std::vector<double> coefB(centres_.size()*trunc_*(r_ + 1), 0.0);

    // loop over data points:
    for(unsigned int i = 0; i < sample.size(); i++)
    {
        // scaled distance between cluster centre and data point:
        double diff = (sample[i] - centres_[idx_[i]])/bw_;
        double expTerm = std::exp(-diff*diff/2.0);

        // power term can be precomputed for efficiency
        // NOTE: this needs double precision to ovoid overflow!
        std::vector<double> powTerm(trunc_ + r_);
        powTerm[0] = 1.0;
        for(unsigned int k = 1; k < trunc_ + r_; k++)
        {
            powTerm[k] = powTerm[k - 1] * diff;
        }
    
        // loop up to truncation number:
        for(unsigned int k = 0; k < trunc_; k++)
        {
            // loop up to derivative order:
            for(unsigned int t = 0; t <= r_; t++)
            {
                coefB[idx_[i]*trunc_*(r_+1) + k*(r_+1) + t] += expTerm
                                                             * powTerm[k + t];
            }
        }
    }

    // precompute the factorial term for efficiency:
    // (need subsequent factorials, avoid repeated calls to factorial())
    // (also compute the inverse here, so that later we can avoid division)
    std::vector<double> facTerm(trunc_);
    facTerm[0] = q_;
    for(unsigned int i = 1; i < trunc_; i++)
    {
        facTerm[i] = facTerm[i - 1]/i;
    }

    // factorial term factored in via separate loop:
    for(unsigned int i = 0; i < centres_.size(); i++)
    {
        for(unsigned int k = 0; k < trunc_; k++)
        {
            for(unsigned int t = 0; t < r_ + 1; t++)
            {
                coefB[i*trunc_*(r_ + 1) + k*(r_ + 1) + t] *= facTerm[k];
            }
        }
    }

    // return coefficient matrix:
    std::vector<real> ret(coefB.begin(), coefB.end());
    return ret;
}


/*!
 * Calculates the scalar coefficient
 *
 * \f[
 *      q = \frac{(-1)^r}{\sqrt{2\pi} n h^{r+1}}
 * \f]
 *
 * which depends on the bandwidth and the number of sample points.
 */
real
GaussianDensityDerivative::setupCoefQ(unsigned int n)
{
    return std::pow(-1, r_) / (std::sqrt(2.0*M_PI)*n*std::pow(bw_, r_ + 1));
}


/*!
 * Determines the cutoff radius according to
 *
 * \f[
 *      R_\text{c} = R_\text{i} + \min(1, 2h\sqrt{\log\big(\sqrt{r!}/\epsilon^\prime\big)})
 * \f]
 *
 * where the scaled error tolerance, \f$ \epsilon^\prime \f$ is determined 
 * using setupScaledTolerance().
 */
real GaussianDensityDerivative::setupCutoffRadius()
{
    // as data is scaled to unit interval, maximum cutoff radius is 1.0:
    real tmp = std::min(1.0, 
               2.0*bw_*std::sqrt(std::log(std::sqrt(rFac_)/epsPrime_)));
    return ri_ + tmp;
}


/*!
 * Computes the appropriately scaled tolerance according to
 *
 * \f[
 *      \epsilon^\prime = \frac{\epsilon}{N*q}
 * \f]
 *
 * where the factor \f$ q \f$ is calculated with setupCoefQ().
 */
real
GaussianDensityDerivative::setupScaledTolerance(unsigned int n)
{
    return eps_/(n*q_);
}


/*
 *
 */
unsigned int
GaussianDensityDerivative::setupTruncationNumber()
{
    // hardcoded limit for truncation number:
    // NOTE: this number will be exponent of distance between point and 
    // cluster, if too large, double will overflow!
    const unsigned int TRUNCMAX = 500;

    // factors constant in the loop:
    double bwSq = bw_*bw_;
    double riSq = ri_*ri_;

    // find lowest truncation number that guarantees error bound:
    for(unsigned int p = 0; p <= TRUNCMAX; p++)
    {
        // calculate error for given truncation number?
        double b = std::min<double>(rc_, 0.5*(ri_ + std::sqrt(riSq + 8.0*p*bwSq)));
        double d = ri_ - b;
        double err = std::sqrt(rFac_)/factorial(p) 
                 * std::pow(ri_*b/bwSq, p) 
                 * std::exp(-0.25*d*d/(bwSq));

        // has error bound been reached?
        if( err < epsPrime_ )
        {
            // plus one for safety:
            return (p + 1);
        }
    }

    // throw exception if error can not be reduced to within desired bound:
    throw std::runtime_error("Could not converge error bound without "
                             "exceeding maximum truncation number p = 500.");
    return -1;
}


/*!
 * Evaluates the hermite polynomial of given order at \f$ x \f$. Uses direct
 * recursion and is potentially not very efficient (but is only used in 
 * reference implementation).
 */
double   
GaussianDensityDerivative::hermite(
        double x, 
        unsigned int r)   
{   
    if( r == 0 )   
    {   
        return 1.0;   
    }   
    else if( r == 1 )   
    {   
        return x;   
    }   
    else   
    {   
        return x*hermite(x, r - 1) - (r - 1)*hermite(x, r - 2);   
    }   
}  


/*!
 * Returns the factorial of the given number. Uses floating point variable to
 * be able to represent factorial of larger numbers beyond the largest integer.
 */
double
GaussianDensityDerivative::factorial(double n)
{
    real f = 1.0;
    for(int i = 1; i <= n; i++)
    {
        f *= i;
    }
    return f;
}


/*!
 * Finds the offset and scaling factor to map both given vectors onto the unit
 * interval.
 */
std::pair<real, real>
GaussianDensityDerivative::getShiftAndScaleParams(
        const std::vector<real> &sample,
        const std::vector<real> &eval)
{
    // sanity checks:
    if( sample.empty() && eval.empty() )
    {
        throw std::runtime_error("Can not calculate shift and scale parameters"
                                 " if both input vectors are empty!");
    }

    // find minimal and maximal values in sample and evaluation vectors:
    auto rangeSample = std::minmax_element(sample.begin(), sample.end());
    auto rangeEval = std::minmax_element(eval.begin(), eval.end());

    // find minimum and maximum over both vectors:
    real minVal = std::min(*rangeSample.first, *rangeEval.first);
    real maxVal = std::max(*rangeEval.second, *rangeEval.second);
    
    // return parameters as pair:
    return std::pair<real, real>(-minVal, 1.0/(maxVal - minVal));
}


/*!
 * Shifts and scales all elements in given vector by a constant offset and
 * factor.
 */
void
GaussianDensityDerivative::shiftAndScale(
        std::vector<real> &vec,
        real shift,
        real scale)
{
    std::for_each(
            vec.begin(), 
            vec.end(), 
            [shift, scale](real &v){v = (v + shift)*scale;}); 
}


/*!
 * Inverts the operation performed by shiftAndScale(). Assuming that the same 
h * shift and scale parameters are used, this will map the input vector back to 
 * its original interval.
 */
void
GaussianDensityDerivative::shiftAndScaleInverse(
        std::vector<real> &vec,
        real shift,
        real scale)
{
    scale = 1.0/scale;
    std::for_each(
            vec.begin(), 
            vec.end(), 
            [shift, scale](real &v){v = v*scale - shift;}); 
}

