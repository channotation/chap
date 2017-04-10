#include <cstdlib>
#include <cstdio>
#include <iostream>

#include <cblas.h>                                                              
#include <lapacke.h>

#include "geometry/cubic_spline_interp.hpp"


/*
 *
 */
CubicSplineInterp::CubicSplineInterp()
    : B_()
{

}


/*
 *
 */
CubicSplineInterp::~CubicSplineInterp()
{

}


/*
 *
 */
int
CubicSplineInterp::interpolate(const std::vector<real> &x, 
                               const std::vector<real> &y)
{
    // sanity check:
    if( x.size() != y.size() )
    {
        std::cout<<"ERROR: x and y vectors must be of same size!"<<std::endl;
        std::abort();
    }

    // dimesnion of linear system:
    int nData = x.size();
    int nSys = nData + 2;

    
    // ASSEMBLE RIGHT HAND SIDE VECTOR
    //-------------------------------------------------------------------------

    // allocate rhs vector:
    real *rhs = new real[nSys];

    // estimate endpoint derivatives:
    // TODO implement this!
    real deriv_lo = 0;
    real deriv_hi = 0;

    // add endpoint conditions to rhs vector:
    rhs[0] = deriv_lo;
    rhs[nSys - 1] = deriv_hi;

    // fill remaining elements of rhs vector:
    for(unsigned int i = 0; i < nData - 1; i++)
    {
        rhs[i + 1] = y[i];
    }


    // ASSEMBLE LEFT HAND SIDE MATRIX
    //-------------------------------------------------------------------------

    // allocate lhs matrix and set elements to zero:
    real *mat = new real[nSys*nSys];

    // prepare knot vector:
    std::vector<real> knotVec;

    // calculate sub diagonal elements:
    std::vector<real> subDiag;
    subDiag.resize(nSys);

    // upper endpoint:
    // TODO
    subDiag[nSys - 2] = 99.3; 

    // internal points:
    for(unsigned int i = 0; i < nSys - 1; i++)
    {
        subDiag[i + 1] = B_(knotVec, degree_, i, x[i]);
    }

    // calculate main diagonal elements:
    std::vector<real> mainDiag;
    mainDiag.resize(nSys);

    // endpoints:
    // TODO
    mainDiag[0] = -99.0;
    mainDiag[nSys - 1] = -99.1;
    
    // internal points:
    for(unsigned int i = 0; i < nData; i++)
    {
        mainDiag[i + 1] = B_(knotVec, degree_, i + 1, x[i]);
    }

    // calculate super diagonal elements:
    std::vector<real> superDiag;
    superDiag.resize(nSys);

    // lower endpoint:
    // TODO
    superDiag[0] = 99.2;

    // internal points:
    for(unsigned int i = 0; i < nData - 1; i++)
    {
        superDiag[i + 1] = B_(knotVec, degree_, i + 2, x[i]);
    }

    // set matrix elements:
    for(unsigned int i = 0; i < nSys; i++)
    {
        for(unsigned int j = 0; j < nSys; i++)
        {
            if( j == i - 1 )
            {
                mat[i*nSys + j] = subDiag[i];
            }
            if( j == i )
            {
                mat[i*nSys + j] = mainDiag[i];
            }
            if( j == i + 1 )
            {
                mat[i*nSys + j] = superDiag[i];
            }
        }
    }
    


    // SOLVE SYSTEM
    //-------------------------------------------------------------------------

    int status = 0;

    // PREPARE OUTPUT
    //-------------------------------------------------------------------------
}



