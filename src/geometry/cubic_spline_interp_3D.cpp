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
#include <stdexcept>
#include <string>

#include <lapacke.h>

#include "geometry/basis_spline.hpp"
#include "geometry/cubic_spline_interp_3D.hpp"


/*!
 * Constructor. 
 */
CubicSplineInterp3D::CubicSplineInterp3D()
{

}


/*!
 * Destructor.
 */
CubicSplineInterp3D::~CubicSplineInterp3D()
{

}


/*!
 * Interpolation interface.
 */
SplineCurve3D
CubicSplineInterp3D::interpolate(std::vector<gmx::RVec> &points,
                                 eSplineInterpBoundaryCondition bc)
{
    // find parameterisation:
    std::vector<real> chordLength = calcChordLength(points);

    // interpolate using this parameterisation:
    return interpolate(chordLength, points, bc);
}


/*!
 * Interpolation interface as convenient operator.
 */
SplineCurve3D
CubicSplineInterp3D::operator()(std::vector<gmx::RVec> &points,
                                eSplineInterpBoundaryCondition bc)
{
    return interpolate(points, bc);
}


/*!
 * Interpolation interface.
 */
SplineCurve3D
CubicSplineInterp3D::interpolate(std::vector<real> &param,
                                 std::vector<gmx::RVec> &points,
                                 eSplineInterpBoundaryCondition bc)
{
    // calculate corresponding knot vector:
    std::vector<real> knotVector = prepareKnotVector(param);

    // get coordinates of each dimension individually:
    std::vector<real> x;
    std::vector<real> y;
    std::vector<real> z;
    for(unsigned int i = 0; i < points.size(); i++)
    {
        x.push_back(points[i][0]);
        y.push_back(points[i][1]);
        z.push_back(points[i][2]);
    }


    // Assmble Left Hand Side Matrix:
    //-------------------------------------------------------------------------

    // dimension of system and number of right hand sides:
    size_t nDat = points.size();
    size_t nSys = nDat + 2;

    // allocate system matrix:
    real subDiag[nSys - 1];
    real mainDiag[nSys];
    real superDiag[nSys - 1];

    // assemble the matrix diagonals:
    assembleDiagonals(knotVector,
                      param,
                      subDiag,
                      mainDiag,
                      superDiag,
                      bc);


    // Assemble Right Hand Side Vector
    //-------------------------------------------------------------------------

    // number of right hand sides is one:
    size_t nRhs = 3;

    // allocate right hand side:
    real rhsMat[nSys * nRhs];
    real rhsX[nSys];
    real rhsY[nSys];
    real rhsZ[nSys];

    // assemble the rhs vectors:
    assembleRhs(param, x, rhsX, bc);
    assembleRhs(param, y, rhsY, bc);
    assembleRhs(param, z, rhsZ, bc);

    // assemble rhs vectors into matrix:
    for(size_t i = 0; i < nSys; i++)
    {
        rhsMat[i] = rhsX[i];
        rhsMat[i + nSys] = rhsY[i];
        rhsMat[i + 2*nSys] = rhsZ[i];
    }

    
    // Solve System
    //-------------------------------------------------------------------------

    // solve tridiagonal system by Gaussian elimination:
    int status = LAPACKE_sgtsv(LAPACK_COL_MAJOR,
                               nSys, 
                               nRhs,
                               subDiag, 
                               mainDiag, 
                               superDiag, 
                               rhsMat,
                               nSys);

    // handle solver failure:
    if( status != 0 )
    {
        throw std::runtime_error("Could not solve tridiagonal system in 3D "
                                 "interpolation. "
                                 "LAPACK error code: "+std::to_string(status));
    } 


    // Prepare Output
    //-------------------------------------------------------------------------

    // create vectorial representation of coefficients:
    std::vector<gmx::RVec> coefs;
    for(size_t i = 0; i < nSys; i++)
    {
        coefs.push_back(gmx::RVec(rhsMat[i],
                                  rhsMat[i + nSys],
                                  rhsMat[i + 2*nSys]));
    }

    // create spline curve object:
    SplineCurve3D SplC(3, knotVector, coefs);

    return SplC;
}


/*!
 * Interpolation interface as convenient operator.
 */
SplineCurve3D
CubicSplineInterp3D::operator()(std::vector<real> &param,
                                std::vector<gmx::RVec> &points,
                                eSplineInterpBoundaryCondition bc)
{
    return interpolate(param, points, bc);
}


/*!
 * Calculates chord length for a given (ordered) point set.
 */
std::vector<real>
CubicSplineInterp3D::calcChordLength(const std::vector<gmx::RVec> &points)
{
    // chord length of first point is zero by definition:
    std::vector<real> chordLength;
    chordLength.push_back(0.0);

    // calculate chord length parameter for all subsequent points:
    for(size_t i = 1; i < points.size(); i++)
    {
        // calculate difference between subsequent points:
        gmx::RVec dVec;
        dVec[XX] = points[i - 1][XX] - points[i][XX];
        dVec[YY] = points[i - 1][YY] - points[i][YY];
        dVec[ZZ] = points[i - 1][ZZ] - points[i][ZZ];
        real diff = std::sqrt(dVec[XX]*dVec[XX] + dVec[YY]*dVec[YY] + dVec[ZZ]*dVec[ZZ]);

        // add to chord length vector:
        chordLength.push_back(chordLength[i-1] + diff);
    }

    return chordLength;
}

