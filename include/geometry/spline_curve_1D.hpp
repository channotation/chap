#ifndef SPLINE_CURVE_1D_HPP
#define SPLINE_CURVE_1D_HPP

#include <utility>
#include <vector>

#include <gtest/gtest_prod.h>   

#include <gromacs/math/vec.h>

#include "geometry/abstract_spline_curve.hpp"


/*!
 * \brief Spline curve in one dimension.
 *
 * This class represents a spline curve in one spatial dimension, i.e. a spline
 * function. In three dimensions, the class SplineCurve3D can be used.
 */
class SplineCurve1D : public AbstractSplineCurve
{
    friend class SplineCurve1DTest;
    FRIEND_TEST(SplineCurve1DTest, SplineCurve1DFindIntervalTest);

    public:
    
        // constructor and destructor:
        SplineCurve1D(
                int degree, 
                std::vector<real> knotVector,
                std::vector<real> ctrlPoints); 
        SplineCurve1D();

        // public interface for curve evaluation:
        real evaluate(
                const real &eval, 
                unsigned int deriv);
        std::vector<real> evaluateMultiple(
                const std::vector<real> &eval, 
                unsigned int deriv);

        // getter function for control points:
        std::vector<real> ctrlPoints() const;

        // compute spline properties:
        real length() const;
        std::pair<real, real> minimum(const std::pair<real, real> &lim);

    private:

        // internal variables:
        std::vector<real> ctrlPoints_;

        // auxiliary functions for evaluation:
        inline real evaluateInternal(const real &eval, unsigned int deriv);
        inline real evaluateExternal(const real &eval, unsigned int deriv);
        inline real computeLinearCombination(const SparseBasis &basis);
};

#endif

