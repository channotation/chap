#include "geometry/basis_spline.hpp"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <cmath>


/*
 * Constructor.
 */
BasisSpline::BasisSpline()
    : evalPoint_(0.0)
    , knotVector_()
{

}


/*
 * Destructor.
 */
BasisSpline::~BasisSpline()
{

}


/*
 * Evaluation function.
 */
inline real 
BasisSpline::evaluate(std::vector<real> &knotVector, 
                      int degree,
                      int interval,
                      real &evalPoint)
{
    // sanity check:
    if( degree < 0 )
    {
        std::cerr<<"ERROR: Polynomial degree must not be negative!"<<std::endl;
        std::abort();
    }

    // initialise state variables:
    evalPoint_ = evalPoint;

    // clear internal knot vector:
    knotVector_.clear();

    // copy internal knots:
    for(int i = 0; i < knotVector.size(); i++)
    {
        knotVector_.push_back(knotVector[i]);
    }

    // handle special case of lying on the upper boundary knot:
    if( evalPoint_ == knotVector_.back() )
    {
        // handle knot multiplicity:
        std::vector<real>::iterator it;
        it = std::lower_bound(knotVector_.begin(), knotVector_.end(), evalPoint);   
        int idx = it - knotVector_.begin();

        // only last basis vector is nozero in this case:
        if( interval == idx - 1 )
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }
    else
    {
        // recursion is handled by separate function:
        return recursion(degree, interval);
    }
}


/*
 * Evaluation function as operator. Refers to evaluate() method.
 */
real
BasisSpline::operator()(std::vector<real> &knotVector,
                        int degree,
                        int interval,
                        real &evalPoint)
{
    // actual calculation is performed by evaluate() method:
    return evaluate(knotVector, degree, interval, evalPoint);
}


/*
 * Function recursively calculates value of basis spline of degree k in 
 * knot interval i, where the evaluation point x and knot vector t are 
 * maintained as state members of the class. This is the Cox-de Boor recursion
 * that is sometimes used as the definition of basis splines.
 */
real
BasisSpline::recursion(int k, int i)
{
    // recursion reaches bottom when polynomial degree is zero:
    if( k == 0 )
    {
        // check if evaluation point lies inside i-th knot interval:
        if( evalPoint_ >= knotVector_[i] && evalPoint_ < knotVector_[i + 1] )
        {
            // evaluation point inside interval:
            return 1.0;
        }
        else
        {
            // evaluation point outside interval:
            return 0.0;
        }
    }
    else
    {
        // calculate numerator and denominator of the two prefactors:
        real frstNum = evalPoint_ - knotVector_[i];
        real frstDen = knotVector_[i + k] - knotVector_[i];
        real scndNum = knotVector_[i + k + 1] - evalPoint_;
        real scndDen = knotVector_[i + k + 1] - knotVector_[i + 1];

        // declare prefactors:
        real frstFac;
        real scndFac;

        // calculate first prefactor:
        if( frstDen <= std::numeric_limits<real>::epsilon() )
        {
            // handle zero division:
            frstFac = 0.0;
        }
        else
        {
            frstFac = frstNum / frstDen;
        }
        
        // calculate second prefactor:
        if( scndDen <= std::numeric_limits<real>::epsilon() )
        {
            // handle zero division:
            scndFac = 0.0;
        }
        else 
        {
            scndFac = scndNum / scndDen;
        }

        // call this function recursively:
        return frstFac*recursion(k - 1, i) + scndFac*recursion(k - 1, i + 1);
    }
}


/*
 * Constructor.
 */
BasisSplineDerivative::BasisSplineDerivative()
{

}


/*
 * Destructor.
 */
BasisSplineDerivative::~BasisSplineDerivative()
{

}


/*
 * Function to evaluate derivative of basis spline of given degree and at given
 * evaluation point. Makes use of the relation
 *
 * dB_[i,k] / dx = k*( B_[i,k-1]/(t_[i+k] - t_i) - B_[i+1, k-1]/(t_[i+k+1] - t_[i+1]) )
 *
 * where t is the knot vector. The value of the lowe degree basis splines is 
 * found using the functor defined above for this purpose.
 */
real
BasisSplineDerivative::evaluate(std::vector<real> &knotVector,
                                int degree,
                                int interval,
                                real &evalPoint)
{
    // create basis spline functor:
    BasisSpline B;
    
    // initialise derivative as zero:
    real deriv = 0.0;    

    // calculate denominators:
    real denomA = (knotVector.at(interval + degree) - knotVector.at(interval)); 
    real denomB = (knotVector.at(interval + degree + 1) - knotVector.at(interval + 1));

    // handle zero division:
    if( std::abs(denomA) >= std::numeric_limits<real>::epsilon() )
    {
        // calculate fist part in derivative sum:
        deriv += B(knotVector, degree - 1, interval, evalPoint) / denomA;
    }
    if( std::abs(denomB) >= std::numeric_limits<real>::epsilon() )
    {
        // calculate second part in derivative sum:
        deriv -= B(knotVector, degree - 1, interval + 1, evalPoint) / denomB;
    }

    // multiply with degree:
    deriv *= degree;

    // return derivative value:
    return deriv;
}


/*
 * Evaluation function as operator. Refers to evaluate() method.
 */
real
BasisSplineDerivative::operator()(std::vector<real> &knotVector,
                                  int degree,
                                  int interval,
                                  real &evalPoint)
{
    // actual calculation is performed by evaluate() method:
    return evaluate(knotVector, degree, interval, evalPoint);
}

