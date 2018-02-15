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


#include "geometry/basis_spline.hpp"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <cmath>


/*!
 * Constructor.
 */
BasisSpline::BasisSpline()
    : evalPoint_(0.0)
    , knotVector_()
{

}


/*!
 * Destructor.
 */
BasisSpline::~BasisSpline()
{

}


/*!
 * This function prepares the evaluation of the basis spline via the 
 * Cox-de-Boor recursion. In particular it handles the edge case of the
 * evaluation point coinciding with the final knot.
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
    // TODO: internal knot vector needed?
    knotVector_.clear();

    // copy internal knots:
    for(size_t i = 0; i < knotVector.size(); i++)
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

        // only last basis vector is nonzero in this case:
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


/*!
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


/*!
 * This function recursively calculates value of basis spline of degree 
 * \f$ k \f$ in knot interval \f$ i \f$, where the evaluation point \f$ x  \f$
 * and knot vector \f$ \mathbf{t} \f$ are maintained as state members of the 
 * class. This implements the Cox-de Boor recursion that is sometimes used as 
 * the definition of basis splines.
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


/*!
 * Destructor.
 */
BasisSplineDerivative::~BasisSplineDerivative()
{

}


/*!
 * Function to evaluate derivative of basis spline of given degree and at given
 * evaluation point. */
real
BasisSplineDerivative::evaluate(std::vector<real> &knotVector,
                                int degree,
                                int interval,
                                real &evalPoint,
                                unsigned int derivOrder)
{
    // create basis spline functor:
    BasisSpline B;
        
    // initialise derivative as zero:
    real deriv = 0.0;   

    // has recursion reached bottom (i.e. first derivative):
    if( derivOrder == 1 )
    {

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
    else
    {
        // calculate denominators::
        real denomA = knotVector.at(interval + degree) - knotVector.at(interval);
        real denomB = knotVector.at(interval + degree + 1) - knotVector.at(interval + 1);

        // handle zero division:
        if( std::abs(denomA) >= std::numeric_limits<real>::epsilon() )
        {
            // calculate first part of derivative sum:
            deriv += evaluate(knotVector, degree - 1, interval, evalPoint, derivOrder - 1) / denomA;
        }
        if( std::abs(denomB) >= std::numeric_limits<real>::epsilon() )
        {
            // calculate second part of derivative sum:
            deriv -= evaluate(knotVector, degree - 1, interval + 1, evalPoint, derivOrder - 1) / denomB;
        }

        // multiply with degree:
        deriv *= (degree);

        // return derivative value:
        return deriv;
    }
}


/*!
 * Evaluation function as operator. Refers to evaluate() method.
 */
real
BasisSplineDerivative::operator()(std::vector<real> &knotVector,
                                  int degree,
                                  int interval,
                                  real &evalPoint,
                                  unsigned int derivOrder)
{
    // actual calculation is performed by evaluate() method:
    return evaluate(knotVector, degree, interval, evalPoint, derivOrder);
}

