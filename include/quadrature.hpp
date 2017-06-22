#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <vector>


/*
 *
 */
class BooleQuadarature
{
    public:

        // quadrature as an operator:
        real operator()(std::vector<real> y);
};

#endif

