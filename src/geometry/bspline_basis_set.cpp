#include "geometry/bspline_basis_set.hpp"


/*
 *
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
    // sanity checks:
    if( deriv > degree )
    {
        // TODO: handle this case more gently, i.e. return zeros as appropriate:
        std::cerr<<"ERROR: deriv > degree is not allowed!"<<std::endl;
        std::abort();
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
    unsigned int nBasis = knots.size() - degree - 1;
    std::vector<real> basisSet(nBasis, 0.0);
    for(size_t i = 0; i < nonzeroBasisElements[deriv].size(); i++)
    {
        basisSet.at(i + knotSpanIdx - degree) = nonzeroBasisElements.at(deriv).at(i);
    }

    // return basis set:
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


/*
 *
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


/*
 *
 */
std::vector<std::vector<real>>
BSplineBasisSet::evaluateNonzeroBasisElements(
        real eval,
        const std::vector<real> &knots,
        unsigned int degree,
        unsigned int deriv,
        unsigned int knotSpanIdx)
{

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

    std::vector<std::vector<real>> a(2, std::vector<real>(degree + 1));

    // loop over function index / basis elements:
    for(size_t i = 0; i <= degree; i++)
    {
        // allocate helper array:
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
/*            std::cout<<"d post first if"<<"  "
                     <<"k = "<<k<<"  "
                     <<"i = "<<i<<"  "
                     <<"d = "<<d<<"  "
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
/*            std::cout<<"d post loop"<<"  "
                     <<"k = "<<k<<"  "
                     <<"i = "<<i<<"  "
                     <<"d = "<<d<<"  "
                     <<std::endl;*/

            // additional summand:
            if( i <= pk )
            {
                a[s2][k] = -a[s1][k-1]/ndu[pk+1][i];
                d += a[s2][k]*ndu[i][pk];
            }
/*            std::cout<<"d post second if"<<"  "
                     <<"k = "<<k<<"  "
                     <<"i = "<<i<<"  "
                     <<"d = "<<d<<"  "
                     <<std::endl;*/

            // assign value of derivative to output matrix:
            ders[k][i] = d;

            std::cout<<"k = "<<k<<"  "
                     <<"i = "<<i<<"  "
                     <<"d = "<<d<<"  "
                     <<"ders[k][i] = "<<ders[k][i]<<"  "
                     <<std::endl;

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

    for(int i = 0; i < ders.size(); i++)
    {
        for(int j = 0; j < ders.at(i).size(); j++)
        {
            std::cout<<ders[i][j]<<"\t";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;

    // return output matrix:
    return ders;
}

