#include "geometry/bspline_basis_element.hpp"


/*
 *
 */
real
BSplineBasisElement::operator()(
        real eval,
        int idx,
        const std::vector<real> &knots,
        unsigned int degree)
{
    // uppermost knot index:
    int m = knots.size() - 1;

    // handle special cases:
    if( idx == 0 && eval == knots[0]  )
    {
        return 1.0;
    }
    if( idx == m - degree - 1 && eval == knots[m] )
    {
        return 1.0;
    }

    // use local support property to simply return zeros:
    if( eval < knots[idx] || eval >= knots[idx + degree + 1] )
    {
        return 0.0;
    }

    // allocate memory for coefficients in triangle table:
    std::vector<real> coefs(degree + 1);

    // initialise degree zero basis functions (bottom of Cox-deBoor recursion):
    for(int i = 0; i <= degree; i++)
    {
        if( eval >= knots[idx + i] &&
            eval < knots[idx + i + 1] )
        {
            coefs[i] = 1.0;
        }
        else
        {
            coefs[i] = 0.0;
        }
    }

    // compute basis spline via triangular table:
    for(int i = 1; i <= degree; i++)
    {
        real saved;
        if( coefs[0] == 0.0 )
        {
            saved = 0.0;
        }
        else
        {
            saved = ((eval - knots[idx])*coefs[0])/(knots[idx+i] - knots[idx]);
        }

        for(int j = 0; j <= degree - i + 1; j++)
        {
            // knots to the left and right:
            real left = knots[idx + j + 1];
            real right = knots[idx + j + i + 1];

            // compute coefficients:
            if( coefs[j + 1] == 0.0 )
            {
                coefs[j] = saved;
                saved = 0.0;
            }
            else
            {
                real tmp = coefs[j + 1]/(right - left);
                coefs[j] = saved + (right - eval)*tmp;
                saved = (eval - left)*tmp;
            }
        }
    }

    // return value of basis function:
    return coefs[0];
}


/*
 *
 */
real
BSplineBasisElement::operator()(
        real eval,
        const std::vector<real> &knots,
        unsigned int degree,
        unsigned int deriv)
{
    return 0;
}







