#include <gtest/gtest.h>


/*!
 * Driver of unit test execution.
 */
int main(int argc, char **argv) {

		// initialise Google testing framework:
		testing::InitGoogleTest(&argc, argv);

		// run all tests:
		return RUN_ALL_TESTS();
}

