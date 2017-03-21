#include <iostream>

#include <gtest/gtest.h>

#include <gromacs/trajectoryanalysis.h>

#include "trajectory-analysis/path_finding_module.hpp"


/*
 * Fixture for unit tests of path finding module.
 */
class PathFindingModuleTest : public ::testing::Test
{
	protected:
		
		// constructor:
		PathFindingModuleTest()
		{
		
		};

};



/*
 *
 */
TEST_F(PathFindingModuleTest, HalloTest)
{
    gmx::RVec initProbePos(0.0, 0.0, 0.0);
    gmx::RVec channelDirectionVector(0.0, 0.0, 1.0);


    gmx::Selection poreSel;

//	PathFindingModule pfm(initProbePos, channelDirectionVector, poreSel);

 //   pfm.findPath();




    ASSERT_TRUE(true);
}

