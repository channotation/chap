#include <gtest/gtest.h>

#include "optim/nelder_mead_module.hpp"


/*
 *
 */
class NelderMeadModuleTest : public ::testing::Test
{

};


/*
 *
 */
TEST_F(NelderMeadModuleTest, NelderMeadModuleRosenbrockTest)
{
    NelderMeadModule nmm;

    std::map<std::string, real> params;
    params["test"] = 1.0;
    nmm.setParams(params);
}

