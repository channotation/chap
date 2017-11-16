#include <algorithm>
#include <cmath>

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
    real scale = setupCoefQ(sample.size());
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
    double d = 0.0;
    for(auto s : sample)
    {
        double diff = (eval - s)/bw_;
        d += std::exp(-0.5*diff*diff)*hermite(diff, r_);
    }
    return d;
}


/*
 *
 */
std::vector<real>
GaussianDensityDerivative::compute_B(std::vector<real> sample)
{
    int r = r_;
    int p = trunc_;
    int K = centres_.size();
    int N = sample.size();
    int num_of_B_terms=centres_.size()*trunc_*(r_+1);  
    real q = q_;
    real h = bw_;
    std::vector<real> px = sample;   
    std::vector<double> B_terms(num_of_B_terms);
    std::vector<real> pClusterCenter = centres_;
    std::vector<unsigned int> pClusterIndex = idx_;

    double *k_factorial;   
    k_factorial=new double[p];   
   
    k_factorial[0]=1;   
    for(int i=1; i<p ;i++){   
        k_factorial[i]=k_factorial[i-1]/i;   
        //printf("%f \n",k_factorial[i]);   
    }   
   
    double *temp3;   
    temp3=new double[p+r];   
   
    for(int n=0; n<K; n++){   
        //printf("Cluster %d ",n);   
        for(int k=0; k<p; k++){   
            for(int m=0; m< r+1; m++){   
                B_terms[(n*p*(r+1))+((r+1)*k)+m]=0.0;;   
                //printf("%f ",B_terms[(n*p*(r+1))+((r+1)*k)+m]);   
            }   
        }   
        //printf("\n");   
    }   
   
    for(int i=0; i<N; i++){   
        int cluster_number=pClusterIndex[i];   
        double temp1=(px[i]-pClusterCenter[cluster_number])/h;   
        double temp2=exp(-temp1*temp1/2);   
        temp3[0]=1;   
        for(int k=1; k<p+r; k++){   
            temp3[k]=temp3[k-1]*temp1;   
        }   
   
        for(int k=0; k<p; k++){   
            for(int m=0; m< r+1; m++){   
                B_terms[(cluster_number*p*(r+1))+((r+1)*k)+m]+=(temp2*temp3[k+m]);   
    

                int tmp = cluster_number*p*(r+1) + k*(r+1) + m;
                if( tmp == 156286 && false )
                {
                    std::cout<<"posterior raykar"<<"  "
                             <<"i = "<<i<<"  "
                             <<"tmp = "<<tmp<<"  "
                             <<"bw = "<<bw_<<"  "
                             <<"eps = "<<eps_<<"  "
                             <<"k = "<<k<<"  "
                             <<"t = "<<m<<"  "
                             <<"e = "<<temp2<<"  "
                             <<"d = "<<temp1<<"  "
                             <<"p = "<<temp3[k+m]<<"  "
                             <<"k! = "<<factorial(k)<<"  "
                             <<"B = "<<B_terms[idx_[i]*trunc_*(r_+1) + (r_+1)*k + m]<<"  "
                             <<"B[-1] = "<<B_terms[idx_[i]*trunc_*(r_+1) + (r_+1)*k + m]<<"  "
                             <<std::endl;
                }
        }   
        }   
    }   
   
    for(int n=0; n<K; n++){   
        //printf("Cluster %d ",n);   
        for(int k=0; k<p; k++){   
            for(int m=0; m< r+1; m++){   
                B_terms[(n*p*(r+1))+((r+1)*k)+m]*=(k_factorial[k]*q);   
                //printf("%f ",B_terms[(n*p*(r+1))+((r+1)*k)+m]);   
            }   
        }   
        //printf("\n");   
    }   
   
   
    delete []k_factorial;   
    delete []temp3;  

    std::vector<real> ret(B_terms.begin(), B_terms.end());
    return ret;
}


/*
 *
 */
std::vector<real>
GaussianDensityDerivative::estimateApprox(
        std::vector<real> &sample,
        std::vector<real> &eval)
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
        deriv.push_back(estimApproxAt(sample, e));
    }
    
    // return vector of density derivative at each evaluation point:
    return deriv;
}


/*
 *
 */
real
GaussianDensityDerivative::estimApproxAt(
        const std::vector<real> & /*sample*/,
        real eval)
{ 
    // upper bound for coefficient loop:
    unsigned int sMax = floor(static_cast<real>(r_)/2.0);

    int skip = 0;
    int count = 0;

    // sum up the terms in approximation of derivative:
    double sumPos = 0.0;
    double sumNeg = 0.0;
    double sum = 0.0;           
    for(unsigned int l = 0; l < centres_.size(); l++)
    {
        // distance from cluster centre:
        double dist = eval - centres_[l];


        count++;
        // skip this iteration if cluster centre beyond cutoff radius:
        if( std::fabs(dist) > rc_ )
        {
            continue; // TODO
            skip++;            
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
                for(unsigned int m = 0; m <= r_ - 2*s; m++)
                {   
                    // add to sum:
                    sum += coefA_[idxA]
                         * coefB_[l*trunc_*(r_ + 1) + (r_ + 1)*k + m]
                         * expTerm
                         * powTerm[k + r_ - 2*s - m];
                  
                    double tmp = coefA_[idxA]
                               * coefB_[l*trunc_*(r_+1) + (r_+1)*k + m]
                               * expTerm
                               * powTerm[k+r_ - 2*s - m];

                    if(std::isnan(tmp))
                    {
                        std::cout<<"bw = "<<bw_<<"  "
                                 <<"eps = "<<eps_<<"  "
                                 <<"trunc = "<<trunc_<<"  "
                                 <<"k = "<<k<<"  "
                                 <<"s = "<<s<<"  "
                                 <<"t = "<<m<<"  "
                                 <<"A = "<<coefA_[idxA]<<"  "
                                 <<"B = "<<coefB_[l*trunc_*(r_+1) + k*(r_+1) + m]<<"  "
                                 <<"d = "<<dist<<"  "
                                 <<"e = "<<expTerm<<"  "
                                 <<"p = "<<powTerm[k+r_ - 2*s - m]<<"  "
                                 <<std::endl;
                    }


                    if(tmp > 0)
                    {
                        sumPos += tmp;
                    }
                    else
                    {
                        sumNeg += tmp;
                    }

                    // increment A-coefficient index:
                    idxA++;
                }   
            }   
        }   
    }   

    double sumAdd = sumNeg + sumPos;
    if( std::fabs(sum - sumAdd) > std::numeric_limits<real>::epsilon() )
    {
        std::cout<<"ATTENTION  "
                 <<"sum = "<<sum<<"  "
                 <<"sumAdd = "<<sumAdd<<"  "
                 <<"diff = "<<std::fabs(sum - sumAdd)<<"  "
                 <<std::endl;
    }
   

    std::cout<<"bw = "<<bw_<<"  "
             <<"eps = "<<eps_<<"  "
             <<"Could have skipped "<<skip<<" out of "<<count<<"."
             <<"   "<<(real)skip/(real)count*100.0<<" %"<<std::endl;

    // return density derivative at evaluation point:
    return sum;
}



/*!
 *
 */
real
GaussianDensityDerivative::estimApproxAtOld2(
        const std::vector<real> &/*sample*/,
        real eval)
{  
    // iteration limit for coefficient loop:
    unsigned int sMax = std::floor(static_cast<real>(r_)/2.0);
       
    // loop over cluster centres:
    real sum = 0.0;
    for(size_t l = 0; l < centres_.size(); l++)
    {  
        // distance from cluster centre:
        real diff = eval - centres_[l];

        // skip this iteration, if cluster centre is beyond cutoff radius:
        if( std::abs(diff) > rc_ )
        {
            continue;
        }

        // scale distance by bandwidth:
        diff /= bw_;
       
        // calculate exponential term:
        real expTerm = exp(-0.5*diff*diff);   

        // the power term can be precomputed for efficiency:
        std::vector<real> powTerm(trunc_ + r_);
        powTerm[0] = 1.0;   
        for(unsigned int i = 1; i < trunc_ + r_; i++)
        {   
            powTerm[i] = powTerm[i - 1] * diff;   
        }   

        // loop up to truncation number:
        for(unsigned int k = 0; k <= trunc_ - 1; k++)
        {   
            // A-coefficient will be accessed in same order as it is created:
            unsigned int idxA = 0;  

            // sums over coefficients:
            for(unsigned int s = 0; s <= sMax; s++)
            {   
                for(unsigned int t = 0; t <= r_ - 2*s; t++)
                {
                    // add to sum:
                    sum += coefA_[idxA] 
                         * coefB_[l*trunc_*(r_+1) + (r_+1)*k + t]
                         * powTerm[k + r_ - 2*s - t]
                         * expTerm;

                    // increment coefficient A index:
                    idxA++;   
                }   
            }   
        }   
    }   

    // return derivative value at this point:
    return sum;
}


/*!
 *
 */
real
GaussianDensityDerivative::estimApproxAtOld(
        const std::vector<real> &/*sample*/,
        real eval)
{
    unsigned int sMax = std::floor(static_cast<real>(r_)/2.0);

    // build sum over centres:
    double sum = 0.0;
    int count = 0;
    for(unsigned int l = 0; l < centres_.size(); l++)
    {
        double d = (eval - centres_[l])/bw_;
        double e = std::exp(-0.5*d*d);
  
        // ignore centres that are more than the cutoff radius from eval point:
        if( std::abs(centres_[l] - eval) > rc_ )
        {
            // FIXME: commenting this back in makes test fail!
            continue;
/*            std::cout<<"l = "<<l<<"  "
                     <<"c = "<<centres_[l]<<"  "
                     <<"eval = "<<eval<<"  "
                     <<"dist = "<<std::abs(centres_[l] - eval)<<"  "
                     <<"rc = "<<rc_<<"  "
                     <<std::endl;*/
        }
/*
        std::vector<real> p(trunc_+r_);
        p[0] = 1.0;
        for(int i = 0; i < trunc_ + r_; i++)
        {
            p[i] = p[i-1]*d;
        }*/

        // sum up to truncation number terms:
        for(unsigned int k = 0; k < trunc_ - 1; k++)
        {
            unsigned int idxA = 0;

            // loops over coefficient matrices:
            for(unsigned int s = 0; s <= sMax; s++)
            {
                for(unsigned int t = 0; t <= r_ - 2*s; t++)
                {
                    real p = std::pow(d, k + r_ - 2*s - t); 
                    real tmp = coefA_.at(idxA)
                         * coefB_.at(l*trunc_*(r_+1) + k*(r_+1) + t) 
                         * e * p;
                    sum += tmp;
/*
                    std::cout<<"l = "<<l<<"  "
                             <<"k = "<<k<<"  "
                             <<"s = "<<s<<"  "
                             <<"t = "<<t<<"  "
                             <<"bw = "<<bw_<<"  "
                             <<"c = "<<centres_[l]<<"  "
                             <<"eval = "<<eval<<"  "
                             <<"d = "<<d<<"  "                           
                             <<"e = "<<e<<"  "
                             <<"p = "<<p<<"  "
                             <<"a = "<<coefA_[idxA]<<"  "
                             <<"b = "<<coefB_[l*trunc_*(r_+1) + k*(r_+1) + t]<<"  "
                             <<"tmp = "<<tmp<<"  "
                             <<std::endl;*/

                    // increment indeces:
                    idxA++;
                    count++;
                }
            }
        }
    }

    std::cout<<"count = "<<count<<std::endl;

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


/*
 *
 */
std::vector<real>
GaussianDensityDerivative::setupCoefA()
{
//    return compute_a();

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


/*
 *
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
    // as data is scaled to unit interval, maximum cutoff radius is 1.0:
    real tmp = std::min(1.0, 
               2.0*bw_*std::sqrt(std::log(std::sqrt(rFac_)/epsPrime_)));
    return ri_ + tmp;
}


/*
 *
 */
real
GaussianDensityDerivative::setupScaledTolerance(unsigned int n) // FIXME remove n?
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
 * shift and scale parameters are used, this will map the input vector back to 
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

