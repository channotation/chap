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
        const std::vector<real> &sample,
        const std::vector<real> &eval)
{
    // allocate output vector of all zeros:
    std::vector<real> deriv;
    deriv.reserve(eval.size());

    // loop over target points:
    for(auto e : eval)
    {
        deriv.push_back(estimApproxAt(sample, e));
    }

    // scale with correct prefactor:
    real scale = std::pow(-1, r_)/(SQRT2PI_*sample.size()*std::pow(bw_, r_+1));
    std::for_each(deriv.begin(), deriv.end(), [scale](real &d){d *= scale;}); 

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
    


    return 0.0;
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
 * Sets derivative order \f$ r>0 \f$.
 */
void
GaussianDensityDerivative::setDerivOrder(unsigned int r)
{
    r_ = r;
}


/*!
 * Sets error bound for approximate method \f$ \epsilon > 0 \f$.
 */
void
GaussianDensityDerivative::setErrorBound(real eps)
{
    eps_ = eps;
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
    unsigned int aNum = sNum*tNum;

    // allocate coeeficient container:
    std::vector<real> coefA;
    coefA.resize(aNum);

    // evaluate the coefficients:
    // TODO: this can be optimised by porecomputing factorials!
    unsigned int idx = 0;
    for(unsigned int s = 0; s <= sMax; s++)
    {        
        for(unsigned int t = 0; t <= tMax; t++)
        {

            std::cout<<"s = "<<s<<"  "<<"sMax = "<<sMax<<"  "
                     <<"t = "<<t<<"  "<<"tMax = "<<tMax<<"  "
                     <<"r = "<<r_<<"  "
                     <<std::endl;
            coefA[idx] = std::pow(-1, s+t) * factorial(r_) 
                        / (std::pow(2, s) * factorial(s) * factorial(t)
                        * factorial(r_ - 2*s - t));
            idx++;
        }
    }

    std::cout<<"ende"<<std::endl;
    return coefA;
}


/*
 *
 */
std::vector<real>
GaussianDensityDerivative::setupCoefARef()
{
    // (inclusive) maximum of indices to coefficient matrix:
    unsigned int sMax = std::floor(r_/2);
    unsigned int tMax = r_;

    unsigned int sNum = sMax + 1;
    unsigned int tNum = tMax + 1;
    unsigned int aNum = sNum*tNum;

    // allocate coeeficient container:
    std::vector<real> coefA;
    coefA.resize(aNum);

    // evaluate the coefficients:
    // TODO: this can be optimised by porecomputing factorials!
    unsigned int idx = 0;
    for(unsigned int s = 0; s <= sMax; s++)
    {        
        for(unsigned int t = 0; t <= tMax; t++)
        {
            std::cout<<"s = "<<s<<"  "
                     <<"t = "<<t<<"  "
                     <<std::endl;
            coefA[idx] = std::pow(-1, s+t) * factorial(r_) 
                        / (std::pow(2, s) * factorial(s) * factorial(t)
                        * factorial(r_ - 2*s - t));
            idx++;
        }
    }


    return coefA;
}


/*!
 * Returns the factorial of the given integer.
 */
unsigned int
GaussianDensityDerivative::factorial(unsigned int n)
{
    std::cout<<"n = "<<n<<std::endl;

    unsigned int f = 1;
    for(int i = 1; i <= n; i++)
    {
        f *= i;
    }
    return f;
}













