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
 *
 */
SplineCurve1D
CubicSplineInterp1D::interpolate(std::vector<real> x,
                                 std::vector<real> f)
{
    
    // generate knot vector:
    std::vector<real> knotVector = prepareKnotVector(x);

    for(int i = 0; i < knotVector.size(); i++)
    {
        std::cout<<"i = "<<i<<"   "
                 <<"knot = "<<knotVector.at(i)<<"  "
                 <<std::endl;
    }

    std::cout<<"blah!"<<std::endl;


    // Assmble Left Hand Side Matrix:
    //-------------------------------------------------------------------------

    // dimension of system and number of right hand sides:
    int nDat = x.size();
    int nSys = nDat + 2;

    // initialise basis spline (derivative) functor:
    BasisSpline B;
    BasisSplineDerivative D;

    // initialse diagonals:
    real subDiag[nSys - 1] = {66, 66, 66, 66, 66, 66, 66};
    real mainDiag[nSys] = {66, 66, 66, 66, 66, 66};
    real superDiag[nSys - 1] = {66, 66, 66, 66, 66, 66};

    for(int i = 0; i < nSys; i++)
    {
        std::cout<<"i = "<<i<<"  "
                 <<"alpha = "<<mainDiag[i]<<"  "
                 <<std::endl;
    }
    std::cout<<"----------------------------------------------------------"<<std::endl;


    // assemble subdiagonal:
    subDiag[nSys - 2] = D(knotVector, degree_, nSys - 2, x.back());
    for(int i = 0; i < nDat; i++)
    {
        subDiag[i] = B(knotVector, degree_, i, x.at(i));
    }

    // assemble main diagonal:
    mainDiag[0] = D(knotVector, degree_, 0, x.front());
    mainDiag[nSys - 1] = D(knotVector, degree_, nSys - 1, x.back());
    for(int i = 0; i < nDat; i++)
    {
        mainDiag[i + 1] = B(knotVector, degree_, i+1, x.at(i));
    }


    // assemble superdiagonal:
    superDiag[0] = D(knotVector, degree_, 1, x.front());
    for(int i = 1; i < nSys - 1; i++)
    {
        superDiag[i] = B(knotVector, degree_, i + 1, x.at(i - 1));
    }



    for(int i = 0; i < nSys; i++)
    {
        std::cout<<"i = "<<i<<"  "
                 <<"alpha = "<<mainDiag[i]<<"  "
                 <<std::endl;
    }
    for(int i = 0; i < nSys - 1; i++)
    {
        std::cout<<"i = "<<i<<"  "
                 <<"beta = "<<subDiag[i]<<"  "
                 <<std::endl;
    }
    for(int i = 0; i < nSys - 1; i++)
    {
        std::cout<<"i = "<<i<<"  "
                 <<"gamma = "<<superDiag[i]<<"  "
                 <<std::endl;
    }



    // Assemble Right Hand Side Vector
    //-------------------------------------------------------------------------

    // number of right hand sides is one:
    int nRhs = 1;
 
    // initialise right hand side:
    real rhsVec[nSys * nRhs];

    // estimate derivatives at endpoints via finite differences:
    rhsVec[0] = ( f.at(1) - f.at(0) ) / ( x.at(1) - x.at(0) );
    rhsVec[nSys - 1] = ( f.at(nDat - 1) - f.at(nDat - 2) ) / ( x.at(nDat - 1) - x.at(nDat - 2) );

    std::cout<<"nSys = "<<nSys<<std::endl;

    // assemble internal points:
    for(int i = 0; i < nDat; i++)
    {
        rhsVec[i + 1] = f.at(i);
    }

     
    
    // Solve System
    //-------------------------------------------------------------------------
   
   
    for(int i  = 0; i < nSys; i++)
    {
        std::cout<<"i = "<<i<<"  "
                 <<"rhs = "<<rhsVec[i]<<std::endl;
    }


    // solve tridiagonal system by Gaussian elimination:

    int status = LAPACKE_sgtsv(LAPACK_COL_MAJOR,
                               nSys, 
                               nRhs,
                               subDiag, 
                               mainDiag, 
                               superDiag, 
                               rhsVec,
                               nSys);

    std::cout<<"status = "<<status<<std::endl;
    // handle solver failure:
    if( status != 0 )
    {
        std::cerr<<"ERROR: Could not solve tridiagonal system!"<<std::endl;
        std::cerr<<"LAPACK error code: "<<status<<std::endl;
        std::abort();
    }

    for(int i  = 0; i < nSys; i++)
    {
        std::cout<<"i = "<<i<<"  "
                 <<"c = "<<rhsVec[i]<<std::endl;
    }


    // Prepare Output
    // -------------------------------------------------------


    for(int i = 0; i < x.size(); i++)
    {
               
        real sum = 0.0;
        for(int j = 0; j < nSys; j++)
        {
            sum += rhsVec[j] * B(knotVector, degree_, j, x.at(i));
        }

        std::cout<<"x = "<<x.at(i)<<"\t"
                 <<"f = "<<f.at(i)<<"\t"
                 <<"s = "<<sum<<std::endl; 
    }




    std::vector<real> ctrlPoints;
    ctrlPoints.resize(nSys);
    for(unsigned int i = 0; i < nSys; i++)
    {
        ctrlPoints.at(i) = rhsVec[i];
    }

    // create spline curve object:
    SplineCurve1D Spl(degree_, knotVector, ctrlPoints);

    // return spline curve:
    return Spl;
}


/*
 *
 */
std::vector<real> 
CubicSplineInterp1D::prepareKnotVector(std::vector<real> &t)
{
    // initialise knot vector:
    std::vector<real> knotVector;

    // add repeats of first knot:
    for(int i = 0; i < degree_; i++)
    {
        knotVector.push_back(t.front());
    }

    // copy support points:
    for(int i = 0; i < t.size(); i++)
    {
        knotVector.push_back(t.at(i));
    }

    // add repeats of last knot:
    for(int i = 0; i < degree_; i++)
    {
        knotVector.push_back(t.back());
    }

    // return knot vector:
    return knotVector;
}

