#include "geometry/linear_spline_interp_1D.hpp"

#include <stdexcept>


/*!
 * Public interpolation interface. Returns a one-dimensional spline curve that
 * linearly interpolates the given set of function values and evaluation 
 * points. Input vectors should be of same length and ordered in the same way.
 */
SplineCurve1D
LinearSplineInterp1D::interpolate(
        const std::vector<real> &x,
        const std::vector<real> &f)
{
    // sanity check:
    if( x.size() != f.size() )
    {
        throw std::logic_error("Input vectors for interpolation must be of "
        "the same size!");
    }

    // degree one for linear interpolation:
    int splineDegree = 1;

    // need to duplicate endpoints to get knot vector:
    std::vector<real> knotVector = x;
    knotVector.push_back(x.back());
    knotVector.insert(knotVector.begin(), x.front());

    // create spline object:
    SplineCurve1D splc(
        splineDegree,
        knotVector,
        f);

    // return the interpolating spline curve:
    return splc;
}


/*!
 * Interface to for linear interpolation as an operator. Simply forwards the
 * arguments to interpolate().
 */
SplineCurve1D
LinearSplineInterp1D::operator()(
        const std::vector<real> &x,
        const std::vector<real> &f)
{
    return interpolate(x, f);
}

