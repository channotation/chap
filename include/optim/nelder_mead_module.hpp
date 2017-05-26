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
        void setObjFun(ObjectiveFunction objFun);
        void setInitGuess(std::vector<real> guess);

        // optimisation and result retrieval:
        void optimise();
        OptimSpacePoint getOptimPoint();
        
    private:

        // control parameters:
        int maxIter_;
        real initShiftFac_;

        // internal parameters:
        real contractionPar_;
        real expansionPar_;
        real reflectionPar_;
        real shrinkagePar_;

        // objective function:
        ObjectiveFunction objFun_;

        // internal optimisation state:
        std::vector<OptimSpacePoint> simplex_;
        OptimSpacePoint centroid_; 

        // comparison functor:
        CompOptimSpacePoints comparison_;

        // internal functions:
        void calcCentroid();
};

#endif

