#ifndef SPLINE_CURVE_1D_HPP
#define SPLINE_CURVE_1D_HPP

#include <vector>

#include <gtest/gtest_prod.h>   

#include <gromacs/math/vec.h>

#include "geometry/abstract_spline_curve.hpp"


/*
 *
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
        ~SplineCurve1D();

        // public interfact for spline evaluation:
        real evaluate(
                real &evalPoint, 
                unsigned int derivOrder, 
                eSplineEvalMethod method);
        real operator()(
                real &evalPoint, 
                unsigned int derivOrder,
                eSplineEvalMethod method);

        // getter function for control points:
        std::vector<real> ctrlPoints() const;

    private:

        // internal variables:
        std::vector<real> ctrlPoints_;
};


#endif

