# Setup Google Test
# -----------------------------------------------------------------------------

# find GTest directory:
find_package(GTest REQUIRED)

# add GTest include directories:
include_directories(${GTEST_INCLUDE_DIRS})



# get list of all source file but ignore main file of chap directory:
file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB_RECURSE TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
list(REMOVE_ITEM SRC_FILES ${PROJECT_SOURCE_DIR}/src/main.cpp)
#message("SRC = " ${SRC_FILES})
#message("TEST_SRC = " ${TEST_SRC_FILES})

# create executable to run all tests:
add_executable(runAllTests runAllTests.cpp ${SRC_FILES} ${TEST_SRC_FILES})
target_link_libraries(runAllTests ${GROMACS_LIBRARIES} ${GTEST_LIBRARIES} ${LAPACK_LIBRARIES} pthread)

# create install target for tets runner:
install(TARGETS runAllTests DESTINATION /sansom/s117/scro2967/repos/chap/bin)




# Move to Test Directories
# -----------------------------------------------------------------------------

#add_subdirectory(statistics-utilities)
#add_subdirectory(trajectoryAnalysis)
