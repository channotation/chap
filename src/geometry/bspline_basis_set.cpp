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

#include "geometry/bspline_basis_set.hpp"


/*!
 * High level public interface for the evaluation of a complete set of B-spline
 * basis functions. Internally this uses evaluateNonZeroBasisElements() to
 * compute the nonzero elements of the basis at the given evaluation point. 
 * The resulting vector of length \f$ p + 1 \f$, where \f$ p \f$ is the spline
 * degree, is then embedded in a vector of length \f$ m - p -1 \f$ where 
 * \f$ m \f$ is the length of the given knot vector, such that each nonzero 
 * basis element is at the correct index. The vector returned by this function
 * then contains the B-spline basis elements \f$ B_{i,p} \f$ with 
 * \f$ i\in[0, m - p - 1] \f$.
 */
SparseBasis
BSplineBasisSet::operator()(
        const real &eval,
        const std::vector<real> &knots,
        unsigned int degree)
{
    // find knot span for evalution point:
    size_t knotSpanIdx = findKnotSpan(eval, knots, degree);

    // calculate the nonzero basis elements:
    std::vector<real> nonzeroBasisElements = evaluateNonzeroBasisElements(
            eval,
            knots,
            degree,
            knotSpanIdx);

    // create sparse basis vector with appropriate indexing:
    SparseBasis basisSet;
    basisSet.reserve(nonzeroBasisElements.size());
    for(size_t i = 0; i < nonzeroBasisElements.size(); i++)
    {
        basisSet[i + knotSpanIdx - degree] = nonzeroBasisElements[i];
    }

    // return nonzero basis functions:
    return basisSet;
}


/*!
 * High level public interface for the evaluation of the derivatives of a 
 * complete set of B-Spline basis functions. Internally this calls 
 * evaluateNonzeroBasisElements() to obtain efficiently evaluate only those 
 * basis elements and derivatives that are not equal to zero. It then selects
 * only the derivatives of the requested degree and embeds them in a vector of
 * zeros such that the nonzero derivatives are located at the correct index.
 * The return value is a vector of the B-spline basis derivative 
 * \f$ B_{i,p}^(n) \f$ with \f$ i\in[0, m - p - 1] \f$, where \f$ p \f$ is the 
 * spline degree and and \f$ n \f$ is the order of the requested derivative 
 * where by convention the zeroth derivative is the basis function itself.
 */
SparseBasis
BSplineBasisSet::operator()(
        real eval,
        const std::vector<real> &knots,
        unsigned int degree,
        unsigned int deriv)
{
    // number of basis elements:
    unsigned int nBasis = knots.size() - degree - 1;

    // reserve memory for sparse basis vector:
    SparseBasis basisSet;
    basisSet.reserve(nBasis);

    // derivative order higher than spline degree:
    if( deriv > degree )
    {
        // simply return vector of all zeros in this case:
        for(unsigned int i = 0; i < nBasis; i++)
        {
            basisSet[i] = 0.0;
        }
        return basisSet;
    }

    // find knot span for evalution point:
    unsigned int knotSpanIdx = findKnotSpan(eval, knots, degree);

    // find nonzero basis elements and their derivatives:
    std::vector<std::vector<real>> nonzeroBasisElements;
    nonzeroBasisElements = evaluateNonzeroBasisElements(
            eval, 
            knots,
            degree,
            deriv,
            knotSpanIdx);

    // pad with zeros to create full length basis vector:
    for(size_t i = 0; i < nonzeroBasisElements[deriv].size(); i++)
    {
        basisSet[i + knotSpanIdx - degree] = nonzeroBasisElements[deriv][i];
    }

    // return basis set:
    return basisSet;
}


/*!
 * Low level evalution of nonzero basis elements. This implements algorithm 
 * A2.2 from The NURBS book and returns a vector of length \f$ p + 1\f$ 
 * containing the nonzero B-spline basis functions \f$ B_{i,p}(x) \f$, where
 * \f$ p \f$ is the spline degree, \f$ x \f$ is the evaluation point, and 
 * \f$ i \in [j-p,j] \f$ is the index of the basis function. The knot span 
 * index \f$ j \f$ can be computed using findKnotSpan().
 */
std::vector<real>
BSplineBasisSet::evaluateNonzeroBasisElements(
        const real &eval,
        const std::vector<real> &knots,
        unsigned int degree,
        unsigned int knotSpanIdx)
{
    // reserve space for nonzero basis functions:
    std::vector<real> nonzeroBasisElements(degree + 1);

    // reserve space for temporary arrays:
    std::vector<real> left(degree + 1);
    std::vector<real> right(degree + 1);

    // calculate all nonzero basis functions:
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
            real tmp = nonzeroBasisElements[j]/(right[j + 1] + left[i - j]);            
            nonzeroBasisElements[j] = saved + right[j + 1]*tmp;
            saved = left[i - j]*tmp;
        }
        nonzeroBasisElements[i] = saved;
    }

    // return vector of nonzero basis elements:
    return nonzeroBasisElements;
}


/*!
 * Finds the index of the knot span for a given evaluation point \f$ x \f$ and
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


/*!
 * Low level evaluation function for nonzero basis elements and nonzero 
 * derivatives. This implements algorithm A2.3 from The NURBS book and returns 
 * a matrix (implemented as vector of vectors) of dimension 
 * \f$ (n+1) \times (p+1) \f$, where the element \f$ (k,i) \f$ contains the 
 * \f$ k \f$-th derivative of the \f$ i \f$-th B-spline basis, i.e. 
 * \f$ B_{i,p}^{(n)}(x) \f$ with \f$ i \in [j-p,j] \f$ and \f$ k \in [0,p]\f$.
 * Note that the knot span index \f$ j \f$ can be computed using findKnotSpan()
 * and the 0-th derivative is by convention the basis function itself.
 *
 * This function does not explicitly check if the condition \f$ n \leq p \f$
 * holds true and this situation should be handled by the calling functions 
 * from the public interface.
 */
std::vector<std::vector<real>>
BSplineBasisSet::evaluateNonzeroBasisElements(
        real eval,
        const std::vector<real> &knots,
        unsigned int degree,
        unsigned int deriv,
        unsigned int knotSpanIdx)
{
    // allocate memory for temporary data matrix:
    std::vector<std::vector<real>> ndu(degree+1, std::vector<real>(degree+1));

    // reserve space for temporary arrays:
    std::vector<real> left;
    left.resize(degree + 1);
    std::vector<real> right;
    right.resize(degree + 1);

    // compute basis functions and keep coefficients required for derivatives:
    ndu[0][0] = 1.0;
    for(size_t i = 1; i <= degree; i++)
    {
        // calculate numerator of left and right terms in recursion formula:
        left[i] = eval - knots[knotSpanIdx + 1 - i];
        right[i] = knots[knotSpanIdx + i] - eval;

        // calculate the actual coefficients:
        real saved = 0.0;
        for(size_t j = 0; j < i; j++)
        {
            ndu[i][j] = right[j + 1] + left[i - j];
            real tmp = ndu[j][i - 1]/ndu[i][j];

            ndu[j][i] = saved + right[j + 1]*tmp;
            saved = left[i - j]*tmp;
        }
        ndu[i][i] = saved;
    }
         

    // allocate matrix of output values:
    std::vector<std::vector<real>> ders(deriv+1, std::vector<real>(degree+1));

    // copy basis functions (zero derivative) into output matrix:
    for(size_t i = 0; i <= degree; i++)
    {
        ders[0][i] = ndu[i][degree];
    }

    // loop over function index / basis elements:
    for(unsigned int i = 0; i <= degree; i++)
    {

        // allocate helper array:
        std::vector<std::vector<real>> a(2, std::vector<real>(degree + 1));
        a[0][0] = 1.0;

        // indices to alternate rows in a:
        int s1 = 0;
        int s2 = 1;

        // loop over derivatives:
        for(unsigned int k = 1; k <= deriv; k++)
        {
            // temporary variable for summing up to derivative at this (k,1);
            real d = 0.0;

            // differences wrt current derivative degree:
            int ik = i - k;
            unsigned int pk = degree - k;

            // if basis element index greater than derivative order, we need to
            // compute a new element in the helper matrix:
            if( i >= k )
            {
                a[s2][0] = a[s1][0]/ndu[pk + 1][ik];
                d = a[s2][0]*ndu[ik][pk];
            }

            // lower index limit for loop over helper array:
            unsigned int rLo;
            if( ik >= -1 )
            {
                rLo = 1;
            }
            else
            {
                rLo = -ik;
            }

            // upper index limit for loop over helper array:
            unsigned int rHi;
            if( i - 1 <= pk )
            {
                rHi = k - 1;
            }
            else
            {
                rHi = degree - i;
            }

            // sum over helper array to compute value of derivative:
            for(unsigned int r = rLo; r <= rHi; r++)
            {
                a[s2][r] = (a[s1][r] - a[s1][r-1])/ndu[pk+1][ik + r];
                d += a[s2][r]*ndu[ik+r][pk];
            }

            // additional summand:
            // NOTE: extra check for 1 != 0 is not in Piegl's algorithm
            if( i <= pk && i != 0 )
            {
                a[s2][k] = -a[s1][k-1]/ndu[pk+1][i];
                d += a[s2][k]*ndu[i][pk];
            }

            // assign value of derivative to output matrix:
            ders[k][i] = d;

            // switch rows:
            int tmp = s1;
            s1 = s2;
            s2 = tmp;
        }
    }

    // multiply derivatives by correct factors (from derivative recursion):
    int fac = degree;
    for(size_t k = 1; k <= deriv; k++)
    {
        for(size_t i = 0; i <= degree; i++)
        {
            ders[k][i] *= fac;
        }
        fac *= (degree - k);
    }

    // return output matrix:
    return ders;
}

