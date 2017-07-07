#ifndef LINEAR_SPLINE_INTERP_1D
#define LINEAR_SPLINE_INTERP_1D

#include <vector>

#include "geometry/spline_curve_1D.hpp"


/*!
 * \brief Functor for performing linear interpolation in one dimension.
 *
 * This class can be used to linearly interpolate a set of function values 
 * \f$ f_i = f(x_i) \f$ defined at a set of evaluation points \f$ x_i \f$. It 
 * will return a spline curve of degree one which goes through all input 
 * points:
 *
 * \f[
 *      s(x_i) = f_i
 * \f]
 *
 * Internally, this is done by setting up a B-spline as
 *
 * \f[
 *      s(x) = \sum_{i=1}^N f_i B_{i,1}(x)
 * \f]
 *
 * where \f$ N \f$ is the number of evaluation points.
 *
 * As a linear spline, the resulting curve will not be smooth, but it will be
 * tight, i.e. it is not prone to overshoot where the curvature of \f$ f(x) \f$
 * is large. This property is useful when probability density  functions, the 
 * value of which may never be negative. It is therefore used in classes 
 * derived from AbstractDensityEstimator.
 */
class LinearSplineInterp1D
{
    public:

        // interpolation interface:
        SplineCurve1D interpolate(
                const std::vector<real> &x,
                const std::vector<real> &f);
        SplineCurve1D operator()(
                const std::vector<real> &x,
                const std::vector<real> &f);

};

#endif

