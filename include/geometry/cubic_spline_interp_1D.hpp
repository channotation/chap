#ifndef CUBIC_SPLINE_INTERP_1D_HPP
#define CUBIC_SPLINE_INTERP_1D_HPP

#include <vector>

#include "geometry/spline_curve_1D.hpp"


enum eSplineInterpBoundaryCondition {eSplineInterpBoundaryHermite, 
                                     eSplineInterpBoundaryNatural};    
enum eSplineInterpEndpoint {eSplineInterpEndpointLo, 
                            eSplineInterpEndpointHi};
enum eDerivEstimate {eDeriv};


/*
 * Functor class for performing cubic spline interpolation in one dimension.
 */
class CubicSplineInterp1D
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

    private:

        // member variables:
        const int degree_ = 3;
        eSplineInterpBoundaryCondition bc_;

        // internal functions:
        std::vector<real> prepareKnotVector(std::vector<real> &x);
        real estimateEndpointDeriv(std::vector<real> &x,
                                   std::vector<real> &f,
                                   eSplineInterpEndpoint endpoint);
};

#endif

