#include <lapacke.h>

#include "geometry/cubic_spline_interp_3D.hpp"


/*
 *
 */
CubicSplineInterp3D::CubicSplineInterp3D()
    : nData_(0)
    , nSys_(0)
{

}


/*
 *
 */
CubicSplineInterp3D::~CubicSplineInterp3D()
{

}


/*
 *
 */
int
CubicSplineInterp3D::interpolate(std::vector<real> t,
                                 std::vector<gmx::RVec> c)
{
    // set system size based on data vector dimension:
    nData_ = c.size();
    nSys_ = nData_ + 2;


    // ESTIMATE DERIVATIVES AT ENDPOINTS
    //-------------------------------------------------------------------------

    // allocate memory for endpoint derivatives:
    real deriv_lo[3];
    real deriv_hi[3];

    // parameter value difference at endpoints:
    real dt_lo = t[1] - t[2];
    real dt_hi = t[nData_] - t[nData_ - 1];
    
    // simple finite difference estimate of derivative:
    deriv_lo[0] = (c[1][0] - c[0][0]) / dt_lo;
    deriv_lo[1] = (c[1][1] - c[0][1]) / dt_lo;
    deriv_lo[2] = (c[1][2] - c[0][2]) / dt_lo;
    deriv_hi[0] = (c[nData_][0] - c[nData_ - 1][0]) / dt_hi;
    deriv_hi[1] = (c[nData_][1] - c[nData_ - 1][1]) / dt_hi;
    deriv_hi[2] = (c[nData_][2] - c[nData_ - 1][2]) / dt_hi;


    // ASSEMBLE RIGHT HAND SIDE MATRIX
    //-------------------------------------------------------------------------
    
    // allocate memory for rhs matrix:
    real *rhs = new real[nSys_*nRhs_];

    // lower endpoint:
    rhs[0*nSys_] = deriv_lo[0];
    rhs[1*nSys_] = deriv_lo[1];
    rhs[2*nSys_] = deriv_lo[2];

    // upper endpoint:
    rhs[1*nSys_ - 1] = deriv_hi[0];
    rhs[2*nSys_ - 1] = deriv_hi[1];
    rhs[3*nSys_ - 1] = deriv_hi[2];

    // loop over internal points:
    for(int i = 1; i < nSys_ - 1; i++)
    {
        rhs[i + 0*nSys_] = c[i][0];
        rhs[i + 1*nSys_] = c[i][1];
        rhs[i + 2*nSys_] = c[i][2];
    }
    
    
    // ASSEMBLE LEFT HAND SIDE MATRIX
    //-------------------------------------------------------------------------

    // allocate memory for diagonals:
    real *subDiag = new real[nSys_ - 1];
    real *mainDiag = new real[nSys_];
    real *superDiag = new real[nSys_ - 1];

    // lower endpoint:
    // TODO!
    mainDiag[0] = -99;
    superDiag[0] = -99;

    // upper endpoint:
    // TODO!
    subDiag[nSys_ - 1] = 99;
    mainDiag[nSys_ - 1] = 99;

    // create knot vector:
    std::vector<real> knotVector = createKnotVector(t);

    // fill in values for intermediate points:
    for(int i = 1; i < nSys_ - 1; i++)
    {
        subDiag[i] = B_(knotVector, degree_, i - 1, t[i - 1]);
        mainDiag[i] = B_(knotVector, degree_, i, t[i - 1]);
        superDiag[i] = B_(knotVector, degree_, i + 1, t[i - 1]);
    }


    // SOLVE LINEAR SYSTEM
    //-------------------------------------------------------------------------

    // solve tridiagonal system by Gaussian elimination:
    int status = LAPACKE_sgtsv(LAPACK_ROW_MAJOR, 
                               nSys_, 
                               nRhs_,
                               subDiag, 
                               mainDiag, 
                               superDiag, 
                               rhs,
                               nSys_);

    // CREATE SPLINE CURVE OBJECT
    //-------------------------------------------------------------------------

    // return status:
    return status;
}


/*
 * Convenience function to generate a valid knot vector from a paramter vector.
 * This is done by simply inserting four copies of the upper and lower 
 * endpoints (i.e. three additional ones).
 */
std::vector<real>
CubicSplineInterp3D::createKnotVector(const std::vector<real> &t)
{
    // create empty knot vector:
    std::vector<real> knotVector;
    
    // append lower endpoint copies:
    knotVector.push_back(t.front());
    knotVector.push_back(t.front());
    knotVector.push_back(t.front());

    // append internal points:
    for(int i = 0; i < t.size(); i++)
    {
        knotVector.push_back(t[i]);
    }

    // append upper endpoint copies:
    knotVector.push_back(t.back());
    knotVector.push_back(t.back());
    knotVector.push_back(t.back());
 
    // return knot vector:
    return knotVector;
}
