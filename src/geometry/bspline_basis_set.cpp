#include "geometry/bspline_basis_set.hpp"


/*
 *
 */
std::vector<real>
BSplineBasisSet::operator()(
        real eval,
        const std::vector<real> &knots,
        unsigned int degree,
        unsigned int deriv)
{
    // find knot span for evalution point:
    size_t knotSpanIdx = findKnotSpan(eval, knots, degree);

    // reserve space for nonzero basis functions:
    std::vector<real> nonzeroBasisElements;
    nonzeroBasisElements.resize(degree + 1);

    // reserve space for temporary arrays:
    std::vector<real> left;
    left.resize(degree + 1);
    std::vector<real> right;
    right.resize(degree + 1);

    // caclualte all nonzero basis functions:
    nonzeroBasisElements[0] = 1.0;
    for(size_t i = 1; i <= degree; i++)
    {
        // calculate numerator of left and right terms in recursion formula:
        left[i] = eval - knots[knotSpanIdx + 1 - i];
        right[i] = knots[knotSpanIdx + i] - eval;

        // calculate value of basis functions:
        real saved = 0.0;
        for(size_t j = 0; j < i; j++)
        {
            real tmp = nonzeroBasisElements[j]/(right.at(j + 1) + left.at(i - j));            
            nonzeroBasisElements[j] = saved + right.at(j + 1)*tmp;
            saved = left.at(i - j)*tmp;
        }
        nonzeroBasisElements[i] = saved;
    }

    // pad with zeros to create full length basis vector:
    unsigned int nBasis = knots.size() - degree - 1;
    std::vector<real> basisSet(nBasis, 0.0);
    for(size_t i = 0; i < nonzeroBasisElements.size(); i++)
    {
        basisSet[i + knotSpanIdx - degree] = nonzeroBasisElements[i];
    }

    // return nonzero basis functions:
    return basisSet;
}


/*!
 * Finds the index of the knot span for a given evalution point \f$ x \f$ and
 * knot vector \f$ \{t\}_{i=1}^m \f$, i.e. finds \f$ i \f$ such that
 * \f$ t_i \leq x < t_{i+1} \f$. The special case of \f$ x = t_m \f$ is handled 
 * as \f$ i = m - p - 1 \f$.
 */
size_t
BSplineBasisSet::findKnotSpan(
        real eval,
        const std::vector<real> &knots,
        unsigned int degree)
{
    // calculate number of basis functions from number of knots and degree:
    unsigned int numBasisFunctions = knots.size() - degree - 1;

    // handle special case of eval point at endpoint:
    if( eval == knots.at(numBasisFunctions)  )
    {
//        std::cout<<"endpoint detected"<<std::endl;
        return numBasisFunctions - 1;
    }

    unsigned int lo = degree;
    unsigned int hi = numBasisFunctions + 1;
    unsigned int mi = (lo + hi)/2;
    while( eval < knots[mi] || eval >= knots[mi + 1] )
    {
        // update interval boundaries:
        if( eval < knots[mi] )
        {
            hi = mi;
        }
        else
        {
            lo = mi;
        }

        // update interval midpoint:
        mi = (lo + hi)/2;
    }

    // return index of knot span:
    return mi;
}


