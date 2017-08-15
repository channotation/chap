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
std::vector<real>
BSplineBasisSet::operator()(
        real eval,
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
std::vector<real>
BSplineBasisSet::operator()(
        real eval,
        const std::vector<real> &knots,
        unsigned int degree,
        unsigned int deriv)
{
    // number of basis elements:
    unsigned int nBasis = knots.size() - degree - 1;

    // derivative order higher than spline degree:
    if( deriv > degree )
    {
        // simply return vector of all zeros in this case:
        return std::vector<real>(nBasis, 0.0);
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
    std::vector<real> basisSet(nBasis, 0.0);
    for(size_t i = 0; i < nonzeroBasisElements[deriv].size(); i++)
    {
        basisSet.at(i + knotSpanIdx - degree) = nonzeroBasisElements.at(deriv).at(i);
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
        real eval,
        const std::vector<real> &knots,
        unsigned int degree,
        unsigned int knotSpanIdx)
{
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

    // return vector of nonzero basis elements:
    return nonzeroBasisElements;
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
 */
std::vector<std::vector<real>>
BSplineBasisSet::evaluateNonzeroBasisElements(
        real eval,
        const std::vector<real> &knots,
        unsigned int degree,
        unsigned int deriv,
        unsigned int knotSpanIdx)
{
//std::cout<<"DERIVATIVE BASED METHOD"<<std::endl;
    //
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

/*
    std::cout<<std::endl;
    std::cout<<"ndu = "<<std::endl;
    for( int i = 0; i <= degree; i++ )
    {
        for( int j = 0; j <= degree; j++ )
        {
            std::cout<<ndu[i][j]<<"\t";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
*/
    // loop over function index / basis elements:
    for(size_t i = 0; i <= degree; i++)
    {

        // allocate helper array:
        std::vector<std::vector<real>> a(2, std::vector<real>(degree + 1));
        a[0][0] = 1.0;

        // indices to alternate rows in a:
        int s1 = 0;
        int s2 = 1;

        // loop over derivatives:
        for(size_t k = 1; k <= deriv; k++)
        {
            // temporary variable for summing up to derivative at this (k,1);
            real d = 0.0;

/*            std::cout<<"d creation"<<"  "
                     <<"k = "<<k<<"  "
                     <<"i = "<<i<<"  "
                     <<"d = "<<d<<"  "
                     <<std::endl;*/

            // differences wrt current derivative degree:
            int ik = i - k;
            int pk = degree - k;

            // if basis element index greater than derivative order, we need to
            // compute a new element in the helper matrix:
            if( i >= k )
            {
                // TODO ERROR: should be = instead of - 
                a[s2][0] = a[s1][0]/ndu[pk + 1][ik];
                d = a[s2][0]*ndu[ik][pk];
            }
/*            std::cout<<"d post frst if"<<"  "
                     <<"k = "<<k<<"  "
                     <<"i = "<<i<<"  "
                     <<"d = "<<d<<"  "
                     <<"idx = "<<knotSpanIdx<<"  "
                     <<std::endl;*/

            // lower index limit for loop over helper array:
            int rLo = 999999;
            if( ik >= -1 )
            {
                rLo = 1;
            }
            else
            {
                rLo = -ik;
            }

            // upper index limit for loop over helper array:
            int rHi = -9999;
            if( i - 1 <= pk )
            {
                rHi = k - 1;
            }
            else
            {
                rHi = degree - i;
            }

            // sum over helper array to compute value of derivative:
            for(int r = rLo; r <= rHi; r++)
            {
                a[s2][r] = (a[s1][r] - a[s1][r-1])/ndu[pk+1][ik + r];
                d += a[s2][r]*ndu[ik+r][pk];

/*
                std::cout<<"  "
                         <<"r = "<<r<<"  "
                         <<"ndu[pk+1][r] = "<<ndu[pk+1][r]<<"  "
                         <<"ndu[ik+r][pk] = "<<ndu[ik+r][pk]<<"  "
                         <<"a[s1][r] = "<<a[s1][r]<<"  "
                         <<"a[s1][r-1] = "<<a[s1][r-1]<<"  "
                         <<"a[s2][r] = "<<a[s2][r]<<"  "
                         <<"d = "<<d<<"  "
                         <<"ik = "<<ik<<"  "
                         <<"pk = "<<pk<<"  "
                         <<std::endl;*/
            }
/*            std::cout<<"d post loop   "<<"  "
                     <<"k = "<<k<<"  "
                     <<"i = "<<i<<"  "
                     <<"d = "<<d<<"  "
                     <<"idx = "<<knotSpanIdx<<"  "
                     <<std::endl;*/

            // additional summand:
            // TODO: there is still an error here for i = 0!
            // 
            // --> simply excluding i = 0 from this if block fixes the first 
            // derivative, but not the second, so I guess that the error is in
            // either a or ndu (either of which should be zero)!
            //
            // --> ndu[0][2] has the same value as in zeroth derivative test
            // and is a basis function which passes the test, so error is 
            // probably in a!
            if( i <= pk & i != 0 )
            {
                a[s2][k] = -a[s1][k-1]/ndu[pk+1][i];
                d += a[s2][k]*ndu[i][pk];
/*                std::cout<<"i = "<<i<<"  "
                         <<"pk = "<<pk<<"  "
                         <<"s1 = "<<s1<<"  "
                         <<"s2 = "<<s2<<"  "
                         <<"a[s1][k-1] = "<<a[s1][k-1]<<"  "
                         <<"ndu[pk+1][i] = "<<ndu[pk+1][i]<<"  "
                         <<"a[s2][k] = "<<a[s2][k]<<"  "
                         <<"ndu[i][pk] = "<<ndu[i][pk]<<"  "
                         <<std::endl;*/
            }
/*            std::cout<<"d post scnd if"<<"  "
                     <<"k = "<<k<<"  "
                     <<"i = "<<i<<"  "
                     <<"d = "<<d<<"  "
                     <<"idx = "<<knotSpanIdx<<"  "
                     <<std::endl;*/

            // assign value of derivative to output matrix:
            ders[k][i] = d;

/*            std::cout<<"k = "<<k<<"  "
                     <<"i = "<<i<<"  "
                     <<"d = "<<d<<"  "
                     <<"ders[k][i] = "<<ders[k][i]<<"  "
                     <<std::endl;*/

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
        // TODO ERROR should be smaller equal
        for(size_t i = 0; i <= degree; i++)
        {
            ders[k][i] *= fac;
        }
        // TODO: bracket this:
        fac *= (degree - k);
    }
/*
    for(int i = 0; i < ders.size(); i++)
    {
        for(int j = 0; j < ders.at(i).size(); j++)
        {
            std::cout<<ders[i][j]<<"\t";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
*/
    // return output matrix:
    return ders;
}

