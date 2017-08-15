#ifndef ABSTRACT_BSPLINE_BASIS_HPP
#define ABSTRACT_BSPLINE_BASIS_HPP

#include <vector>

#include <gtest/gtest.h>

#include <gromacs/utility/real.h>


/*
 *
 */
class AbstractBSplineBasis
{
    public:

        // method for finding the correct knot span:
        size_t findKnotSpan(
                real eval,
                const std::vector<real> &knots,
                unsigned int degree);

};

#endif

