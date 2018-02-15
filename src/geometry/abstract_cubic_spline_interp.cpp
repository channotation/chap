// CHAP - The Channel Annotation Package
// 
// Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
// Stephen J. Tucker
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


#include <iostream>

#include "geometry/abstract_cubic_spline_interp.hpp"


/*!
 * Constructor.
 */
AbstractCubicSplineInterp::AbstractCubicSplineInterp()
{

}


/*!
 * Destructor.
 */
AbstractCubicSplineInterp::~AbstractCubicSplineInterp()
{

}


/*!
 * This function assembles the nonzero entries of the system matrix occurring
 * in spline interpolation. Currently, only Hermite boundary conditions are 
 * implemented.
 */
void
AbstractCubicSplineInterp::assembleDiagonals(std::vector<real> &knotVector,
                                             std::vector<real> &x,
                                             real *subDiag,
                                             real *mainDiag,
                                             real *superDiag,
                                             eSplineInterpBoundaryCondition bc)
{
    // dimension of system:
    int nDat = x.size();
    int nSys = nDat + 2;

    // initialise basis spline (derivative) functor:
    BasisSpline B;
    BasisSplineDerivative D;

    // handle boundary conditions:
    if( bc == eSplineInterpBoundaryHermite )
    {
        int firstOrderDeriv = 1;
        real xLo = x.front();
        real xHi = x.back();

        // lower boundary:
        mainDiag[0] = D(knotVector, degree_, 0, xLo, firstOrderDeriv);
        superDiag[0] = D(knotVector, degree_, 1, xLo, firstOrderDeriv);
 
        // higher boundary:
        mainDiag[nSys - 1] = D(knotVector, degree_, nSys - 1, xHi, firstOrderDeriv);
        subDiag[nSys - 2] = D(knotVector, degree_, nSys - 2, xHi, firstOrderDeriv);
    }
    else if( bc == eSplineInterpBoundaryNatural )
    {
        std::cerr<<"ERROR: Only Hermite boundary conditions implemented!"<<std::endl;
        std::abort();
    }
    else
    {
        std::cerr<<"ERROR: Only Hermite boundary conditions implemented!"<<std::endl;
        std::abort();
    }

    // assemble subdiagonal:
    for(int i = 0; i < nDat; i++)
    {
        subDiag[i] = B(knotVector, degree_, i, x.at(i));
    }

    // assemble main diagonal:
    for(int i = 0; i < nDat; i++)
    {
        mainDiag[i + 1] = B(knotVector, degree_, i+1, x.at(i));
    }

    // assemble superdiagonal:
    for(int i = 1; i < nSys - 1; i++)
    {
        superDiag[i] = B(knotVector, degree_, i + 1, x.at(i - 1));
    }
}


/*!
 * This function assembles the right hand side vector of the tridiagonal 
 * system occurring in cubic spline interpolation. Currently only Hermite 
 * boundary conditions are implemented.
 */
void
AbstractCubicSplineInterp::assembleRhs(std::vector<real> &x,
                                       std::vector<real> &f,
                                       real *rhsVec,
                                       eSplineInterpBoundaryCondition bc)
{
    // get system size:
    int nDat = x.size();
    int nSys = nDat + 2;

    // handle bondary conditions:
    if( bc == eSplineInterpBoundaryHermite )
    {
        // lower boundary:
        rhsVec[0] = estimateEndpointDeriv(x, 
                                          f,
                                          eSplineInterpEndpointLo,
                                          eSplineInterpDerivParabolic);

        // higher boundary:
        rhsVec[nSys - 1] = estimateEndpointDeriv(x, 
                                                 f, 
                                                 eSplineInterpEndpointHi,
                                                 eSplineInterpDerivParabolic);
    }
    else if( bc == eSplineInterpBoundaryNatural )
    {
        // TODO: implement this
        std::cerr<<"ERROR: Only Hermite boundary conditions implemented!"<<std::endl;
        std::abort();
 
        // lower boundary:
        rhsVec[0] = 0.0;

        // higher boundary:
        rhsVec[nSys - 1] = 0.0; 
    }   
    else
    {
        std::cerr<<"ERROR: Only Hermite boundary conditions are implemented!"<<std::endl;
        std::abort();
    }

    // assemble internal points:
    for(int i = 0; i < nDat; i++)
    {
        rhsVec[i + 1] = f.at(i);
    }
}


/*!
 * Internal helper function for creating a knot vector from a vector of input 
 * data. The knot vector is essentially a copy of the data vector with its 
 * first and last element each repeated four times.
 */
std::vector<real> 
AbstractCubicSplineInterp::prepareKnotVector(std::vector<real> &x)
{
    // initialise knot vector:
    std::vector<real> knotVector;

    // add repeats of first knot:
    for(int i = 0; i < degree_; i++)
    {
        knotVector.push_back(x.front());
    }

    // copy support points:
    for(size_t i = 0; i < x.size(); i++)
    {
        knotVector.push_back(x.at(i));
    }

    // add repeats of last knot:
    for(int i = 0; i < degree_; i++)
    {
        knotVector.push_back(x.back());
    }

    // return knot vector:
    return knotVector;
}


/*!
 * This function estimates the derivatives of \f$ f(x) \f$ at the endpoints of
 * the given data range. Estimation can be done with a simple finite difference
 * or via a parabolic fit with a ghost node assumption.
 */
real
AbstractCubicSplineInterp::estimateEndpointDeriv(std::vector<real> &x,
                                                 std::vector<real> &f,
                                                 eSplineInterpEndpoint endpoint,
                                                 eSplineInterpDerivEstimate method)
{
    // get overall number of data points:
    unsigned int nDat = x.size();

    // which finite difference should be used?
    if( method == eSplineInterpDerivParabolic )
    {
        // declare helper variables for local points:
        real xDeltaLo;
        real xDeltaHi;
        real fDeltaLo;
        real fDeltaHi;

        // estimate derivative at lower or higher endpoint?
        if( endpoint == eSplineInterpEndpointLo )
        {
            xDeltaLo = x.at(0) - x.at(2);
            xDeltaHi = x.at(1) - x.at(0);
            fDeltaLo = (f.at(0) - f.at(2))/xDeltaLo;
            fDeltaHi = (f.at(1) - f.at(0))/xDeltaHi;
        }
        else // no else if to avoid un-initialised variable error
        {
            xDeltaLo = x.at(nDat - 1) - x.at(nDat - 2); 
            xDeltaHi = x.at(nDat - 3) - x.at(nDat - 1);
            fDeltaLo = (f.at(nDat - 1) - f.at(nDat - 2))/xDeltaLo;
            fDeltaHi = (f.at(nDat - 3) - f.at(nDat - 1))/xDeltaHi;
        }

        // parabolic estimate of endpoint derivative:
        return (xDeltaLo*fDeltaHi + xDeltaHi*fDeltaLo)/(xDeltaLo + xDeltaHi);
    }
    else 
    {
        // declare helper variables for local points:
        real fHi;
        real fLo;
        real xHi;
        real xLo;

        // estimate derivative at lower or upper endpoint?
        if( endpoint == eSplineInterpEndpointLo )
        {  
           fHi = f.at(1);
           xHi = x.at(1);
           fLo = f.at(0);
           xLo = f.at(0);
        }
        else // no else if to avoid uninitialised variable error
        {
           fHi = f.at(nDat - 1);
           xHi = x.at(nDat - 1);
           fLo = f.at(nDat - 2);
           xLo = f.at(nDat - 2);        
        }

        // simple estimate of endpoint derivative:
        return (fHi - fLo)/(xHi - xLo);
    }
}

