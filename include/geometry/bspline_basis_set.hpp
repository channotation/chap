#ifndef BSPLINE_BASIS_SET_HPP
#define BSPLINE_BASIS_SET_HPP

#include <vector>

#include <gtest/gtest.h>

#include <gromacs/utility/real.h> 


/*
 *
 */
class BSplineBasisSet
{
    friend class BSplineBasisSetTest;
    FRIEND_TEST(BSplineBasisSetTest, BSplineBasisSetKnotSpanTest);

    public:

        // public interface for evaluation:
        std::vector<real> operator()(
                real eval, 
                const std::vector<real> &knots, 
                unsigned int degree);
        std::vector<std::vector<real>> operator()(
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
        /*
        inline std::vector<real> evaluateNonzeroBasisElements(
                real eval,
                const std::vector<real> &knots,
                unsigned int degree,
                unsigned int deriv);*/
};

#endif

