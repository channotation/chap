#ifndef CUBIC_SPLINE_INTERP_3D_HPP
#define CUBIC_SPLINE_INTERP_3D_HPP

#include <vector>

#include "geometry/abstract_cubic_spline_interp.hpp"
#include "geometry/spline_curve_3D.hpp"


/*
 *
 */
class CubicSplineInterp3D : public AbstractCubicSplineInterp
{
    public:

        // constructor and destructor:
        CubicSplineInterp3D();
        ~CubicSplineInterp3D();

        // interpolation interface:
        SplineCurve3D interpolate(const std::vector<gmx::RVec> &points,
                                  eSplineInterpBoundaryCondition bc);
        SplineCurve3D operator()(const std::vector<gmx::RVec> &points,
                                 eSplineInterpBoundaryCondition bc);

    private:

        // internal helpers:
        std::vector<real> calcChordLength(const std::vector<gmx::RVec> &points);
};

#endif

