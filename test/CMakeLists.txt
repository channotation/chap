# Setup Google Test
# -----------------------------------------------------------------------------

# find GTest directory:
find_package(GTest REQUIRED)

# add GTest include directories:
include_directories(${GTEST_INCLUDE_DIRS})


# Move to Test Directories
# -----------------------------------------------------------------------------

add_subdirectory(trajectoryAnalysis)
