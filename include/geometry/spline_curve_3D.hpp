#ifndef SPLINE_CURVE_3D_HPP
#define SPLINE_CURVE_3D_HPP

#include <vector>

#include <gtest/gtest_prod.h>   

#include <gromacs/utility/real.h> 
#include <gromacs/math/vec.h>

#include "geometry/abstract_spline_curve.hpp"


/*!
 * \brief Spline curve in three dimensions.
 *
 * This class represents spline curves in three dimensions. It internally 
 * maintains a set of knots and three-dimensional control points, from which 
 * the expression
 *
 * \f[
 *      \mathbf{S}(s) = \sum_i \mathbf{c}_i B_{i, p}(s)
 * \f]
 *
 * can be used to evaluate the curve and its derivative. This in turn also
 * gives access to the curve's differential properties such as its length(),
 * tangentVec(), and normalVec(). 
 *
 * The method arcLengthParam() can be used to change the internal
 * representation of the curve such that it is parameterised by arc length. 
 */
class SplineCurve3D : public AbstractSplineCurve
{
    public:
      
        // constructor and destructor:
        SplineCurve3D(int degree,
                      std::vector<real> knotVector,
                      std::vector<gmx::RVec> ctrlPoints);
        SplineCurve3D();

        // public interface for curve evaluation:
        gmx::RVec evaluate(const real &eval, unsigned int deriv);

        // reparameterisation methods:
        void arcLengthParam();
        // map points onto curve:
        double pointSqDist(gmx::RVec point, double eval);
        gmx::RVec cartesianToCurvilinear(const gmx::RVec &cartPoint);

        // calculate differential properties of curve:
        real length(const real &lo, const real &hi);
        real length();
        gmx::RVec tangentVec(const real &eval);
        gmx::RVec normalVec(const real/* &evalPoints*/);        
        real speed(const real &eval);

        // utilities for accessing arc length at the control points:
        std::vector<real> ctrlPointArcLength();
        real frstPointArcLength();
        real lastPointArcLength();
 
        // getter functions:
        std::vector<gmx::RVec> ctrlPoints() const;
        
    private:

        // internal variables:
        std::vector<gmx::RVec> ctrlPoints_;
        std::vector<gmx::RVec> refPoints_;

        // arc length lookup table utilities:
        bool arcLengthTableAvailable_;
        std::vector<real> arcLengthTable_;

        // curve evaluation utilities:
        inline gmx::RVec evaluateInternal(const real &eval, unsigned int deriv);
        inline gmx::RVec evaluateExternal(const real &eval, unsigned int deriv);
        inline gmx::RVec computeLinearCombination(const SparseBasis &basis);

        // curve length utilities:
        inline real arcLengthBoole(const real &lo, const real &hi);
        void prepareArcLengthTable();
        
        // arc length reparameterisation utilities:
        inline real arcLengthToParam(real &arcLength);
        inline bool arcLengthToParamTerm(real lo, real hi, real tol);
        inline real arcLengthToParamObj(real lo, real hi, real target);

        //
        unsigned int closestSplinePoint(const gmx::RVec &point);
        gmx::RVec projectionInInterval(
                const gmx::RVec &point,
                const real &lo,
                const real &hi);
        gmx::RVec projectionInExtrapRange(
                const gmx::RVec &point,
                const real &ds);
};

#endif

