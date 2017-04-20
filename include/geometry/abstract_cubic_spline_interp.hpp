#ifndef ABSTRACT_CUBIC_SPLINE_INTERP_HPP
#define ABSTRACT_CUBIC_SPLINE_INTERP_HPP

#include <vector>

#include <gromacs/math/vec.h>

#include "geometry/basis_spline.hpp"


enum eSplineInterpBoundaryCondition {eSplineInterpBoundaryHermite, 
                                     eSplineInterpBoundaryNatural};    
enum eSplineInterpEndpoint {eSplineInterpEndpointLo, 
                            eSplineInterpEndpointHi};
enum eSplineInterpDerivEstimate {eSplineInterpDerivSimple,
                                 eSplineInterpDerivParabolic};


/*
 *
 */
class AbstractCubicSplineInterp
{
    public:

        // constructor and destructor:
        AbstractCubicSplineInterp();
        ~AbstractCubicSplineInterp();

    protected:

        // member variables:
        const int degree_ = 3;
        eSplineInterpBoundaryCondition bc_;

        // internal helper functions:
        void assembleDiagonals(std::vector<real> &knotVector,
                               std::vector<real> &x,
                               real *subDiag,
                               real *mainDiag,
                               real *superDiag,
                               eSplineInterpBoundaryCondition bc);
        void assembleRhs(std::vector<real> &x,
                         std::vector<real> &f,
                         real *rhsVec,
                         eSplineInterpBoundaryCondition bc);
        std::vector<real> prepareKnotVector(std::vector<real> &x);
        real estimateEndpointDeriv(std::vector<real> &x,
                                   std::vector<real> &f,
                                   eSplineInterpEndpoint endpoint,
                                   eSplineInterpDerivEstimate method);
};

#endif

