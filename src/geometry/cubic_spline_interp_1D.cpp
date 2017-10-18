#include <iostream>

#include <lapacke.h>

#include "geometry/basis_spline.hpp"
#include "geometry/cubic_spline_interp_1D.hpp"


/*
 * Constructor.
 */
CubicSplineInterp1D::CubicSplineInterp1D()
{

}



/*
 * Destructor.
 */
CubicSplineInterp1D::~CubicSplineInterp1D()
{

}


/*
 * Public interface for interpolation. Function takes a one-dimensional data 
 * cloud of (x_i, f(x_i)) points (each in their separate vectors) and finds the
 * cubic spline curve that interpolates between these data points such that:
 *
 *      s(x_i) = f(x_i)
 *
 * Currently only Hermite endpoint conditions are implemented. The relevant 
 * linear system is solved via Gaussian elimination (using LAPACKE) and the 
 * result is returned as a splien curve object.
 *
 * TODO: Implement natural boundary conditions.
 */
SplineCurve1D
CubicSplineInterp1D::interpolate(std::vector<real> &x,
                                 std::vector<real> &f,
                                 eSplineInterpBoundaryCondition bc)
{
    std::cout<<"cubic spline interp BEGIN"<<std::endl;

    // sanity check:
    if( x.size() != f.size() )
    {
        std::cerr<<"ERROR: x and f vectors must be of same size!"<<std::endl;
        std::abort();
    }

    // set boundary condition:
    bc_ = bc;

    // generate knot vector:
    std::vector<real> knotVector = prepareKnotVector(x);

    // Assmble Left Hand Side Matrix:
    //-------------------------------------------------------------------------

    // dimension of system and number of right hand sides:
    size_t nDat = x.size();
    size_t nSys = nDat + 2;

    // allocate memory for lhs matrix diagonal entries:
    real subDiag[nSys - 1];
    real mainDiag[nSys];
    real superDiag[nSys - 1];

    // assemble the matrix diagonals:
    assembleDiagonals(knotVector,
                      x,
                      subDiag,
                      mainDiag,
                      superDiag,
                      bc);


    // Assemble Right Hand Side Vector
    //-------------------------------------------------------------------------

    // number of right hand sides is one:
    size_t nRhs = 1;
 
    // initialise right hand side:
    real rhsVec[nSys * nRhs];

    // assemble the rhs vector:
    assembleRhs(x, f, rhsVec, bc);

    
    // Solve System
    //-------------------------------------------------------------------------
  
    // FIXME there are NaNs here that prevent solving the system!
    /*
    std::cout<<"RHS = "<<std::endl;
    for(int i = 0; i < nSys*nRhs; i++)
    {
        std::cout<<rhsVec[i]<<std::endl;
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
                               rhsVec,
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

    // create vector of control points:
    std::vector<real> ctrlPoints;
    ctrlPoints.resize(nSys);
    for(size_t i = 0; i < nSys; i++)
    {
        ctrlPoints.at(i) = rhsVec[i];
    }

    // create spline curve object:
    SplineCurve1D Spl(degree_, knotVector, ctrlPoints);

    // return spline curve:
    return Spl;
}


/*
 * Interpolation interface conveniently defined as operator.
 */
SplineCurve1D
CubicSplineInterp1D::operator()(std::vector<real> &x,
                                std::vector<real> &f,
                                eSplineInterpBoundaryCondition bc)
{
    // actual computation is handled by interpolate() method:
    return interpolate(x, f, bc);
}

