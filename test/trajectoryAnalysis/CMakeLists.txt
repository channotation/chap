# get all *.cpp files recursively
file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/src/trajectoryAnalysis/*.cpp)
file(GLOB_RECURSE TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/trajectoryAnalysis/*.cpp)

# create executable to run all tests:
add_executable(runAllTests runAllTests.cpp ${SRC_FILES} ${TEST_SRC_FILES})
target_link_libraries(runAllTests ${GROMACS_LIBRARIES} ${GTEST_LIBRARIES} pthread)

# create install target for tets runner:
install(TARGETS runAllTests DESTINATION /sansom/s117/scro2967/repos/chap/bin)
