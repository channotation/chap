#include <algorithm>
#include <iostream>

#include "geometry/abstract_spline_curve.hpp"




/*
 *
 */
int
AbstractSplineCurve::findInterval(real evalPoint)
{
    // find upper and lower bounds:
    std::vector<real>::iterator it_hi = std::upper_bound(knotVector_.begin(), 
                                                         knotVector_.end(), 
                                                         evalPoint);

    std::vector<real>::iterator it_lo = std::lower_bound(knotVector_.begin(), 
                                                         knotVector_.end(), 
                                                         evalPoint);
    // get corresponding indeces:
    int idx_hi = it_hi - knotVector_.begin() - 1;
    int idx_lo = it_lo - knotVector_.begin() - 1;



    return idx_hi;
}


/*
 *
 */
real
AbstractSplineCurve::deBoorRecursion(int r, 
                                     int i, 
                                     real &evalPoint,
                                     const std::vector<real> &ctrlCoefs)
{
    // has bottom of recursion been reached?
    if( r == 0 )
    {
        // return control point:
        return ctrlCoefs[i];
    }
    else
    {
        // calculate prefactor denominator:
        real alpha = knotVector_[i + degree_ + 1 - r];
        
        // handle zero division case:
        if( alpha < std::numeric_limits<real>::epsilon() )
        {
            // define 0/0 = 0:
            alpha = 0.0;
        }
        else
        {
            // calculate de Boor coefficient:
            alpha = (evalPoint - knotVector_[i]) / alpha;
        }

        // call recursion:
        return (1.0 - alpha) * deBoorRecursion(r - 1, i - 1, evalPoint, ctrlCoefs) 
                     + alpha * deBoorRecursion(r - 1, i, evalPoint, ctrlCoefs);
    }
}


