#include <cmath>
#include <iostream>

#include <gromacs/math/3dtransforms.h>

#include "path-finding/mock_pore_maker.hpp"


/*
 *
 */
std::vector<gmx::RVec>
MockPoreMaker::makePore(real poreLength,
                        real poreCentreRadius, 
                        real poreVdwRadius,
                        gmx::RVec poreCentre,
                        int alongAxis)
{
    // calculate angle required for overlapping vdW spheres:
    real phi = std::acos(1.0 - std::pow(poreVdwRadius, 2.0)/2.0/std::pow(poreCentreRadius, 2.0));
    int nStepsAround = std::ceil(2.0*PI_/phi);

    // step length and number of steps along the length of the pore:
    real stepLengthAlong = poreVdwRadius/2.0;
    int nStepsAlong = std::ceil(poreLength/stepLengthAlong) + 1;
   
    // place particles on cylinder surface centered around origin:
    std::vector<gmx::RVec> particleCentres;
    for(int i = 0; i < nStepsAlong; i++)
    {
        for(int j = 0; j < nStepsAround; j++)
        {
            // add new particle to pore:
            gmx::RVec particleCentre;
            particleCentre[XX] = poreCentreRadius*std::cos(phi*j);
            particleCentre[YY] = poreCentreRadius*std::sin(phi*j);
            particleCentre[ZZ] = i*stepLengthAlong - poreLength/2.0;
            particleCentres.push_back(particleCentre);
        }
    }

    // rotate in desired way:
    mat4 rotMat;
    vec4 v;
    if( alongAxis == XX )
    {
        // rotate around y-axis:
        gmx_mat4_init_rotation(YY, PI_/2.0, rotMat);

    }
    else if( alongAxis == YY )
    {   
        // rotate around x-axis:
        gmx_mat4_init_rotation(YY, PI_/2.0, rotMat);
    }
    else if( alongAxis = ZZ )
    {
        // nothing to do:
        gmx_mat4_init_rotation(ZZ, 0.0, rotMat);
    }
    else
    {
        std::cerr<<"ERROR: AlongAxis must be one of XX, YY, or ZZ!"<<std::endl;
        std::abort();
    }
    

    for(unsigned int i = 0; i < particleCentres.size(); i++)
    {
        gmx_mat4_transform_point(rotMat, particleCentres[i], v);
        particleCentres[i][XX] = v[XX];
        particleCentres[i][YY] = v[YY];
        particleCentres[i][ZZ] = v[ZZ];
    }



    // shift particles away from origin:
    for(unsigned int i = 0; i < particleCentres.size(); i++)
    {
        particleCentres[i][XX] += poreCentre[XX];
        particleCentres[i][YY] += poreCentre[YY];
        particleCentres[i][ZZ] += poreCentre[ZZ];
    }

    // return mock pore particle positions:
    return particleCentres;
}

