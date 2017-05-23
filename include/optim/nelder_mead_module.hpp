#ifndef NELDER_MEAD_MODULE_HPP
#define NELDER_MEAD_MODULE_HPP

#include <map>
#include <string>

#include <gtest/gtest.h>

#include <gromacs/utility/real.h>

#include "optim/optimisation.hpp"


/*
 *
 */
class NelderMeadModule
{
    friend class NelderMeadModuleTest;
    FRIEND_TEST(NelderMeadModuleTest, NelderMeadModuleRosenbrockTest);

    public:

        // constructor and destructor:
        NelderMeadModule();
        ~NelderMeadModule();

        // setting parameters and initial point:
        void setParams(std::map<std::string, real> params);

        // optimisation and result retrieval:
        void optimise();
        
    private:

        // parameters:
        int maxIter_;


        // internal optimisation state:
        std::vector<OptimSpacePoint> crntSimplex_;
};

#endif

