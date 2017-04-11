#ifndef BASIS_SPLINE_HPP
#define BASIS_SPLINE_HPP

#include <vector>

#include <gromacs/utility/real.h> 

/*
 * Functor class to evaluate basis spline.
 */
class BasisSpline
{
    public:
        
        // constructor and destructor:
        BasisSpline();
        ~BasisSpline();

        // evalution function and operator:
        real evaluate(std::vector<real> &knotVector, 
                      int degree,
                      int interval,
                      real &evalPoint);
        real operator()(std::vector<real> &knotVector,
                        int degree,
                        int interval,
                        real &evalPoint);

    private:

        // state variables:
        real evalPoint_;
        std::vector<real> knotVector_;

        // bspline recursion:
        real recursion(int k, int i);
};


/*
 * Functor class to evaluate basis spline derivative.
 */
class BasisSplineDerivative
{
    public:

        // constructor and destructor:
        BasisSplineDerivative();
        ~BasisSplineDerivative();

        // evaluation function and operator:
        real evaluate(std::vector<real> &knotVector,
                      int degree,
                      int interval,
                      real &evalPoint);
        real operator()(std::vector<real> &knotVector,
                        int degree,
                        int interval, 
                        real &evalPoint);
};

#endif

