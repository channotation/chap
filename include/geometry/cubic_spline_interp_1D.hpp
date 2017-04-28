#ifndef CUBIC_SPLINE_INTERP_1D_HPP
#define CUBIC_SPLINE_INTERP_1D_HPP

#include <vector>

#include "geometry/abstract_cubic_spline_interp.hpp"
#include "geometry/spline_curve_1D.hpp"


/*!
 * \brief Functor class for performing cubic spline interpolation in one 
 * dimension.
 */
class CubicSplineInterp1D : public AbstractCubicSplineInterp
{
    public:
        
        // constructor and destructor:
        CubicSplineInterp1D();
        ~CubicSplineInterp1D();

        // interpolation interface:
        SplineCurve1D interpolate(std::vector<real> &x,
                                  std::vector<real> &f,
                                  eSplineInterpBoundaryCondition bc);
        SplineCurve1D operator()(std::vector<real> &x,
                                 std::vector<real> &f,
                                 eSplineInterpBoundaryCondition bc);
};

#endif

