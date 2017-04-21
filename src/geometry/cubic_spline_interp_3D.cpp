#include <iostream>

#include <lapacke.h>

#include "geometry/basis_spline.hpp"
#include "geometry/cubic_spline_interp_3D.hpp"


/*
 * Constructor. 
 */
CubicSplineInterp3D::CubicSplineInterp3D()
{

}


/*
 * Destructor.
 */
CubicSplineInterp3D::~CubicSplineInterp3D()
{

}


/*
 *
 */
SplineCurve3D
CubicSplineInterp3D::interpolate(const std::vector<gmx::RVec> &points,
                                 eSplineInterpBoundaryCondition bc)
{
    // find parameterisation:
    std::vector<real> chordLength = calcChordLength(points);

    // calculate corresponding knot vector:
    std::vector<real> knotVector = prepareKnotVector(chordLength);

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
    int nDat = points.size();
    int nSys = nDat + 2;

    // allocate system matrix:
    real subDiag[nSys - 1];
    real mainDiag[nSys];
    real superDiag[nSys - 1];

    // assemble the matrix diagonals:
    assembleDiagonals(knotVector,
                      chordLength,
                      subDiag,
                      mainDiag,
                      superDiag,
                      bc);


    // Assemble Right Hand Side Vector
    //-------------------------------------------------------------------------

    // number of right hand sides is one:
    int nRhs = 3;

    // allocate right hand side:
    real rhsMat[nSys * nRhs];
    real rhsX[nSys];
    real rhsY[nSys];
    real rhsZ[nSys];

    // assemble the rhs vectors:
    assembleRhs(chordLength, x, rhsX, bc);
    assembleRhs(chordLength, y, rhsY, bc);
    assembleRhs(chordLength, z, rhsZ, bc);

    // assemble rhs vectors into matrix:
    for(unsigned int i = 0; i < nSys; i++)
    {
        rhsMat[i] = rhsX[i];
        rhsMat[i + nSys] = rhsY[i];
        rhsMat[i + 2*nSys] = rhsZ[i];
    }

    
    // Solve System
    //-------------------------------------------------------------------------
/* 
    std::cout<<"chord length"<<std::endl;
    for(int i = 0; i < chordLength.size(); i++)
    {
        std::cout<<chordLength[i]<<std::endl;
    }
    std::cout<<std::endl;
 

    std::cout<<"knot vector"<<std::endl;
    for(int i = 0; i < knotVector.size(); i++)
    {
        std::cout<<knotVector[i]<<std::endl;
    }
    std::cout<<std::endl;



    std::cout<<"RHS"<<std::endl;
    for(int i = 0; i < nSys; i++)
    {
        std::cout<<rhsX[i]<<"  "
                 <<rhsY[i]<<"  "
                 <<rhsZ[i]<<std::endl;
    }
    std::cout<<std::endl;


    std::cout<<"subDiag"<<std::endl;
    for(int i = 0; i < nSys - 1; i++)
    {
        std::cout<<subDiag[i]<<std::endl;
    }
    std::cout<<std::endl;

    std::cout<<"mainDiag"<<std::endl;
    for(int i = 0; i < nSys; i++)
    {
        std::cout<<mainDiag[i]<<std::endl;
    }
    std::cout<<std::endl;

    std::cout<<"subperDiag"<<std::endl;
    for(int i = 0; i < nSys - 1; i++)
    {
        std::cout<<superDiag[i]<<std::endl;
    }
    std::cout<<std::endl;
*/
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
        std::cerr<<"ERROR: Could not solve tridiagonal system!"<<std::endl;
        std::cerr<<"LAPACK error code: "<<status<<std::endl;
        std::abort();
    } 


    // Prepare Output
    //-------------------------------------------------------------------------

    // create vectorial representation of coefficients:
    std::vector<gmx::RVec> coefs;
    for(unsigned int i = 0; i < nSys; i++)
    {
        coefs.push_back(gmx::RVec(rhsMat[i],
                                  rhsMat[i + nSys],
                                  rhsMat[i + 2*nSys]));
    }

    // create spline curve object:
    SplineCurve3D SplC(3, knotVector, coefs);

    return SplC;
}


/*
 * Interpolation interface as convenient operator.
 */
SplineCurve3D
CubicSplineInterp3D::operator()(const std::vector<gmx::RVec> &points,
                                eSplineInterpBoundaryCondition bc)
{
    return interpolate(points, bc);
}


/*
 * Calculates chord length for a given (ordered) point set.
 */
std::vector<real>
CubicSplineInterp3D::calcChordLength(const std::vector<gmx::RVec> &points)
{
    // chord length of first point is zero by definition:
    std::vector<real> chordLength;
    chordLength.push_back(0.0);

    // calculate chord length parameter for all subsequent points:
    for(unsigned int i = 1; i < points.size(); i++)
    {
        // calculate difference between subsequent points:
        gmx::RVec dVec;
        dVec[0] = points[i - 1][0] - points[i][0];
        dVec[1] = points[i - 1][1] - points[i][1];
        dVec[2] = points[i - 1][2] - points[i][2];
        real diff = std::sqrt(dVec[0]*dVec[0] + dVec[1]*dVec[1] + dVec[2]*dVec[2]);

        // add to chord length vector:
        chordLength.push_back(chordLength[i-1] + diff);
    }

    return chordLength;
}

