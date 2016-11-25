# get all *.cpp files recursively
file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/src/statistics-utilities/*.cpp)
file(GLOB_RECURSE TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/statistics-utilities/*.cpp)

# create executable to run all tests:
add_executable(runStatsTests runStatsTests.cpp ${SRC_FILES} ${TEST_SRC_FILES})
target_link_libraries(runStatsTests ${GROMACS_LIBRARIES} ${GTEST_LIBRARIES} ${LAPACK_LIBRARIES} pthread)

# create install target for tets runner:
install(TARGETS runStatsTests DESTINATION /sansom/s117/scro2967/repos/chap/bin)
