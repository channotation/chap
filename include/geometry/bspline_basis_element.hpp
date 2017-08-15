#ifndef BSPLINE_BASIS_ELEMENT_HPP
#define BSPLINE_BASIS_ELEMENT_HPP

#include <vector>

#include <gtest/gtest.h>

#include <gromacs/utility/real.h>

#include "geometry/abstract_bspline_basis.hpp"


/*
 *
 */
class BSplineBasisElement : public AbstractBSplineBasis
{
    public:

        //
        real operator()(
                real eval,
                int idx,
                const std::vector<real> &knots,
                unsigned int degree);

        //
        real operator()(
                real eval,
                const std::vector<real> &knots,
                unsigned int degree,
                unsigned int deriv);
};

#endif

