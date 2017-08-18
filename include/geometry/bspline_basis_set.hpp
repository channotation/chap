#ifndef BSPLINE_BASIS_SET_HPP
#define BSPLINE_BASIS_SET_HPP

#include <type_traits>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include <gromacs/utility/real.h> 


/*!
 * Shorthand notation for a sparse basis vector given as map of indeces to 
 * values.
 */
typedef std::unordered_map<unsigned int, real> SparseBasis;



/*!
 * \brief Functor class for evaluating complete set of B-spline basis
 * functions and corresponding derivatives.
 *
 * This class implements algorithms for efficiently evaluating a complete 
 * B-spline basis set. This is done by using algorithms A2.2 and A2.3 from The
 * NURBS Book by Piegl and Tiller. This algorithms eveluate each required
 * coefficient only once and are hence significantly more efficient than
 * evaluating each basis function with an individual Cox-deBoor recursion. It
 * is usually the most sensible algorithm to use if e.g. a spline curve needs
 * to be evaluated.
 *
 * In cases where only a single basis function (or derivative thereof) is 
 * required, BSplineBasisElement may be more efficient.
 */
class BSplineBasisSet
{
    friend class BSplineBasisSetTest;
    FRIEND_TEST(BSplineBasisSetTest, BSplineBasisSetKnotSpanTest);

    public:

        // public interface for evaluation:
        SparseBasis operator()(
                const real &eval, 
                const std::vector<real> &knots, 
                unsigned int degree);
   
        // public interface for evaluation with derivatives:
        SparseBasis operator()(
                real eval, 
                const std::vector<real> &knots, 
                unsigned int degree, 
                unsigned int deriv);

    private:

        // method for finding the correct knot span:
        size_t findKnotSpan(
                real eval,
                const std::vector<real> &knots,
                unsigned int degree);

        // method for evaluating the nonzero elements of basis:
        inline std::vector<real> evaluateNonzeroBasisElements(
                const real &eval,
                const std::vector<real> &knots,
                unsigned int degree,
                unsigned int knotSpanIdx);

        // method for evaluating nonzero elements of basis (derivatives):
        inline std::vector<std::vector<real>> evaluateNonzeroBasisElements(
                real eval,
                const std::vector<real> &knots,
                unsigned int degree,
                unsigned int deriv,
                unsigned int knotSpanIdx);
};

#endif

