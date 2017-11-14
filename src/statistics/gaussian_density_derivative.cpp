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
    real d = 0.0;
    for(auto s : sample)
    {
        real diff = (eval - s)/bw_;
        d += std::exp(-0.5*diff*diff)*hermite(diff, r_);
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

//    auto ss = getShiftAndScaleParams(sample, eval);
//    shiftAndScale(sample, ss.first, ss.second);
//    shiftAndScale(eval, ss.first, ss.second);
/*
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
    bw_ *= scale;*/

    // calculate space partitioning (this is data dependent, bc bw_ is scaled):
    centres_ = setupClusterCentres();
    idx_ = setupClusterIndices(sample);

    // compute data dependent coefficients:
    q_ = setupCoefQ(sample.size());
    epsPrime_ = setupScaledTolerance(sample.size());
    rc_ = setupCutoffRadius();
    trunc_ = setupTruncationNumber();
    coefB_ = setupCoefB(sample);

    std::cout<<"trunc = "<<trunc_<<std::endl;

    // loop over target points:
    std::vector<real> deriv;
    deriv.reserve(eval.size());
    for(auto e : eval)
    {
        deriv.push_back(estimApproxAt(sample, e));
    }

        
    // scale data and evaluation points back to original interval:
//    shiftAndScaleInverse(sample, ss.first, ss.second);
//    shiftAndScaleInverse(eval, ss.first, ss.second);

    // scale back data and evaluation point:
    /*
    scale = 1.0/scale;
    std::for_each(sample.begin(), sample.end(), [scale](real &s){ s *= scale; });
    std::for_each(eval.begin(), eval.end(), [scale](real &e){ e *= scale; });
    std::for_each(sample.begin(), sample.end(), [shift](real &s){s += shift;});
    std::for_each(eval.begin(), eval.end(), [shift](real &e){e += shift; });
    bw_ *= scale;*/

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

    // build sum over centres:
    double sum = 0.0;
    for(unsigned int l = 0; l < centres_.size(); l++)
    {
        // ignore centres that are more than the cutoff radius from eval point:
        if( std::abs(centres_[l] - eval) > rc_ )
        {
            continue;
        }

        // sum up to truncation number terms:
        for(unsigned int k = 0; k < trunc_ - 1; k++)
        {
            unsigned int idxA = 0;

            // loops over coefficient matrices:
            for(unsigned int s = 0; s <= sMax; s++)
            {
                for(unsigned int t = 0; t <= r_ - 2*s; t++)
                {
                    real d = (eval - centres_[l])/bw_;
                    real e = std::exp(-0.5*d*d);
                    real p = std::pow(d, k + r_ - 2*s - t); 
                    sum += coefA_.at(idxA)
                         * coefB_.at(l*trunc_*(r_+1) + k*(r_+1) + t) 
                         * e * p;

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
        const std::vector<real> &sample)
{
    // allocate coefficient matrix as NaN:
    std::vector<real> coefB(centres_.size()*trunc_*(r_+1), 0.0);

    // loop over data points:
    for(unsigned int i = 0; i < sample.size(); i++)
    {
        // scaled distance between cluster centre and data point:
        real d = (sample[i] - centres_[idx_[i]])/bw_;
        real e = std::exp(-0.5*d*d);

        
        std::vector<real> p(trunc_ + r_);
        p[0] = 1.0;
        for(unsigned int k = 1; k < trunc_ + r_; k++)
        {
            p[k] = p[k-1]*d;
        }
    

        // loop up to truncation number:
        for(int k = 0; k < trunc_; k++)
        {
            // loop up to derivative order:
            for(int t = 0; t <= r_; t++)
            {
                coefB[idx_[i]*trunc_*(r_+1) + k*(r_+1) + t] += e*std::pow(d, k+t)/factorial(k);
//                coefB[idx_[i]*trunc_*(r_+1) + k*(r_+1) + t] += e*p.at(k+t)/factorial(k);
            }
        }
    }

    // scale all coefficients by common prefactor:
    real fac = q_;
    std::for_each(coefB.begin(), coefB.end(), [fac](real &b){b*=fac;});

    // return coefficient matrix:
    return coefB;
}


/*
 *
 */
std::vector<real>
GaussianDensityDerivative::compute_B(const std::vector<real> &sample)
{
	int p = trunc_;
	int K = centres_.size();
	int r = r_;
    int num_of_B_terms = K*p*(r+1);
	double h = bw_;
	double q = q_;
	std::vector<unsigned int> pClusterIndex = idx_;
	std::vector<real> pClusterCenter = centres_; 
	int N = sample.size();
	std::vector<real> px = sample;
 
    //printf("K=%d p=%d r=%d num_of_B_terms=%d\n",K,p,r,num_of_B_terms);   
   
    double *B_terms=new double[num_of_B_terms];   
   
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

	std::vector<real> coefB;
	for(int i = 0; i < num_of_B_terms; i++)
	{
		coefB.push_back(B_terms[i]);
	} 

	return coefB;
}


/*
 *
 */
std::vector<real>
GaussianDensityDerivative::Evaluate(
        std::vector<real> eval,
        std::vector<real> sample)   
{   
   
    //int num_of_influential_neighbors=(int)ceil(ry/rx);   
    //printf("Num of influential right or left neighbors = %d\n",num_of_influential_neighbors);   
  
    int r = r_;
    int p = trunc_;
    int M = eval.size();
    int K = centres_.size();

    real h = bw_;
    real two_h_square = 2*h*h;

    real rx = 1.0/K;
    real ry = rc_; 

    std::vector<real> py = eval;
    std::vector<real> px = sample;
    std::vector<real> pD(eval.size());
    std::vector<real> pClusterCenter = centres_;
    std::vector<real> a_terms = coefA_;
    std::vector<real> B_terms = coefB_;



    double *temp3;   
    temp3=new double[p+r];   
   
   
    for(int j=0; j<M; j++){   
        pD[j]=0.0;   
           
        int target_cluster_number=std::min((int)floor(py[j]/rx),K-1);   
        double temp1=py[j]-pClusterCenter[target_cluster_number];   
        double dist=abs(temp1);   
        while (dist <= ry && target_cluster_number <K && target_cluster_number >=0){   
   
            //printf("j=%d y=%f Influential cluster=%d\n",j,py[j],target_cluster_number);   
            //Do something   
            double temp2=exp(-temp1*temp1/two_h_square);   
            double temp1h=temp1/h;   
            temp3[0]=1;   
            for(int i=1; i<p+r; i++){   
                temp3[i]=temp3[i-1]*temp1h;   
            }   
   
            for(int k=0; k<=p-1; k++){   
                int dummy=0;   
                for(int l=0; l <= (int)floor((double)r/2); l++){   
                    for(int m=0; m <= r-(2*l); m++){   
                        pD[j]=pD[j]+(a_terms[dummy]*B_terms[(target_cluster_number*p*(r+1))+((r+1)*k)+m]*temp2*temp3[k+r-(2*l)-m]);   
                        dummy++;   
                    }   
                }   
            }   
            //   
                   
   
            target_cluster_number++;   
            temp1=py[j]-pClusterCenter[target_cluster_number];   
            dist=abs(temp1);   
        }   
   
        target_cluster_number=std::min((int)floor(py[j]/rx),K-1)-1;   
        if (target_cluster_number >=0){   
            double temp1=py[j]-pClusterCenter[target_cluster_number];   
            double dist=abs(temp1);   
            while (dist <= ry && target_cluster_number <K && target_cluster_number >=0){   
                //printf("j=%d y=%f Influential cluster=%d\n",j,py[j],target_cluster_number);   
                //Do something   
                double temp2=exp(-temp1*temp1/two_h_square);   
                double temp1h=temp1/h;   
                temp3[0]=1;   
                for(int i=1; i<p+r; i++){   
                    temp3[i]=temp3[i-1]*temp1h;   
                }   
   
                for(int k=0; k<=p-1; k++){   
                    int dummy=0;   
                    for(int l=0; l <= (int)floor((double)r/2); l++){   
                        for(int m=0; m <= r-(2*l); m++){   
                            pD[j]=pD[j]+(a_terms[dummy]*B_terms[(target_cluster_number*p*(r+1))+((r+1)*k)+m]*temp2*temp3[k+r-(2*l)-m]);   
                            dummy++;   
                        }   
                    }   
                }   
                //   
                target_cluster_number--;   
                temp1=py[j]-pClusterCenter[target_cluster_number];   
                dist=abs(temp1);   
            }   
        }   
   
    }   
   
    delete []temp3;


    return pD;
} 


/*
 *
 */
std::pair<std::vector<real>, std::vector<int>>
GaussianDensityDerivative::space_sub_division(std::vector<real> sample)
{   
    int K = centres_.size();
    int N = sample.size();
    real rx = ri_;
    std::vector<real> px = sample;
    std::vector<int> pClusterIndex(sample.size());
    std::vector<real> pClusterCenter(K);


    // 1. Cluster Centers   
   
//    pClusterCenter=new double[K];   
    for(int i=0; i<K; i++){   
        pClusterCenter[i]=(i*rx)+(rx/2.0);   
        //printf("%f\n",pClusterCenter[i]);

    }   
   
    //2. Allocate each source to the corresponding interval   
   
//    pClusterIndex=new int[N];   
    for(int i=0; i<N; i++){   
        pClusterIndex[i]=std::min((int)floor(px[i]/rx),K-1);   
        //printf("x=%f Cluster=%d\n",px[i],pClusterIndex[i]);  
    }   
   
    
    std::pair<std::vector<real>, std::vector<int>> ret(pClusterCenter, pClusterIndex);
//    ret.first = pClusterCenter;
//    ret.second = pClusterIndex;

    return ret;

}   


/*
 *
 */
std::vector<real>
GaussianDensityDerivative::compute_a()
{   
    int r = r_;

    double r_factorial=(double)factorial(r);   
    //printf("%f \n",r_factorial);   
   
    double *l_constant;   
    l_constant=new double[((int)floor((double)r/2))+1];   
    l_constant[0]=1;   
    for(int l=1; l <= (int)floor((double)r/2); l++){   
        l_constant[l]=l_constant[l-1]*(-1.0/(2*l));   
        //printf("%f \n",l_constant[l]);   
    }   
   
    double *m_constant;   
    m_constant=new double[r+1];   
    m_constant[0]=1;   
    for(int m=1; m <= r; m++){   
        m_constant[m]=m_constant[m-1]*(-1.0/m);   
        //printf("%f \n",m_constant[m]);   
    }   
   
    int num_of_a_terms=0;   
    for(int l=0; l <= (int)floor((double)r/2); l++){   
        for(int m=0; m <= r-(2*l); m++){            
            num_of_a_terms++;   
        }   
    }   
   
    //printf("r=%d num_of_a_terms=%d\n",r,num_of_a_terms);   
   
    std::vector<real> a_terms(num_of_a_terms);   
    int k=0;   
    for(int l=0; l <= (int)floor((double)r/2); l++){   
        for(int m=0; m <= r-(2*l); m++){   
            a_terms[k]=(l_constant[l]*m_constant[m]*r_factorial)/((double)factorial(r-(2*l)-m));   
            //printf("%f \n",a_terms[k]);   
            k++;               
        }   
    }   
    delete []l_constant;   
    delete []m_constant;   
  
    return a_terms;
} 





/*
 *
 */
std::vector<real>
GaussianDensityDerivative::EvaluateDirect(std::vector<real> sample, std::vector<real> eval)   
{   
    int M = eval.size();
    int N = sample.size();
    int r = r_;
    real h = bw_;

    std::vector<real> pD(M);
    std::vector<real> py = eval;
    std::vector<real> px = sample;

    double two_h_square=2*h*h;   
    double pi=3.14159265358979;   
    double q=(pow(-1,r))/(sqrt(2*pi)*N*(pow(h,(r+1))));   
   
    for(int j=0; j<M; j++)   
    {   
        pD[j]=0.0;   
   
        for(int i=0; i<N; i++)   
        {   
            double temp=py[j]-px[i];   
            double norm=temp*temp;   
               
            pD[j] = pD[j]+(hermite(temp/h,r)*exp(-norm/two_h_square));             
   
        }   
        pD[j]=pD[j]*q;   
    }   


    return pD;
}   
   



/*
 *
 */
void
GaussianDensityDerivative::choose_parameters(std::vector<real> sample,
	std::vector<real> eval)
{   
    real h = bw_;
    real h_square = h*h;
    real two_h_square = 2*h_square;
    real r = r_;
    real eps = eps_;

    const int P_UL = 100;

    // 1. rx --> interval length.   
   
    real rx=h/2;   
   
    // 2. K  --> number of intervals.   
   
    int K=(int)ceil(1.0/rx);   
    rx=1.0/K;   
    double rx_square=rx*rx;   
   
    // 3. rr --> cutoff radius.   
   
    //double r_term=pow(2.0,r/2)*sqrt((double)factorial(r));   
    double r_term=sqrt((double)factorial(r));   
 
    const real R = 1.0;
    real rr=std::min<real>(R,2*h*sqrt(log(r_term/eps))); // FIXME why min here?  
   
    // 4. ry --> cluster cutoff radius.   
    double ry=rx+rr;   
   
    // 5. p  --> truncation number.   
   
    int p=0;   
    double error=1;   
    double temp=1;   
    double comp_eps=eps/r_term;   
           
    while((error > comp_eps) & (p <= P_UL)){   
        p++;   
        double b=std::min(((rx+sqrt((rx_square)+(8*p*h_square)))/2),ry);   
        double c=rx-b;   
        temp=temp*(((rx*b)/h_square)/p);   
        error=temp*(exp(-(c*c)/2*two_h_square));        
    }      
    p=p+1;    
   
    //printf("h=%f r=%d eps=%f K=%d rx=%f rr=%f ry=%f p=%d\n",h,r,eps,K,rx,rr,ry,p);   
   
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
    return bw_/2.0 + tmp;
}


/*
 *
 */
real
GaussianDensityDerivative::setupScaledTolerance(unsigned int n)
{
    return eps_/std::sqrt(factorial(r_));    
}


/*
 *
 */
unsigned int
GaussianDensityDerivative::setupTruncationNumber()
{
    // hardcoded limit for truncation number:
    const unsigned int TRUNCMAX = 500;

    // factors constant in the loop:
    real bwSq = bw_*bw_;
    real riSq = ri_*ri_;

    // find lowest truncation number that guarantees error bound:
    for(int p = 0; p <= TRUNCMAX; p++)
    {
        // calculate error for given truncation number?
        real b = std::min<real>(rc_, ri_ + 0.5*std::sqrt(riSq + 8.0*p*bwSq));
        real d = ri_ - b;
        real err = std::sqrt(rFac_)/factorial(p) 
                 * std::pow(ri_*b/bwSq, p) 
                 * std::exp(-0.5*d*d/(bwSq));

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
real   
GaussianDensityDerivative::hermite(
        real x, 
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
real
GaussianDensityDerivative::factorial(real n)
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

