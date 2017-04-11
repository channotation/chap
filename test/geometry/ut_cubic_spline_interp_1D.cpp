#include <gtest/gtest.h>

#include "geometry/cubic_spline_interp_1D.hpp"



/*
 *
 */
class CubicSplineInterp1DTest : public ::testing::Test
{


};



/*
 *
 */
TEST_F(CubicSplineInterp1DTest, CubicSplineInterp1DFirstTest)
{
    // 
    std::vector<real> x = {-2.0, -1.0, -0.5,  0.0, 0.5,  1.0, 2.0};
    std::vector<real> f = { 8.0,  1.0,  0.25, 0.0, 0.25, 1.0, 8.0};


    CubicSplineInterp1D Interp;

    SplineCurve1D Spl = Interp.interpolate(x, f);
   

    for(unsigned int i = 0; i < x.size(); i++)
    {
        real val = Spl.evaluate(x.at(i), eSplineEvalDeBoor);
        std::cout<<"x = "<<x.at(i)<<"  "
                 <<"s(x) = "<<val<<"  "
                 <<std::endl;
        ASSERT_NEAR(f.at(i), val, std::numeric_limits<real>::epsilon());
    }


}
