# Setup Google Test
# -----------------------------------------------------------------------------

# build tests only if GTest is available:
find_package(GTest)
if( GTEST_FOUND )

    # create target for running make check:
    add_custom_target(check ${CMAKE_CTEST_COMMAND} -V)

    # include Google test directories:
    include_directories(${GTEST_INCLUDE_DIRS})

    # get list of all source files but ignore main file of chap directory:
    file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
    file(GLOB_RECURSE TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
    list(REMOVE_ITEM SRC_FILES ${PROJECT_SOURCE_DIR}/src/main.cpp)

    # add executable to run all tests and link libraries:
    add_executable(runAllTests ${TEST_SRC_FILES} ${SRC_FILES})
    target_link_libraries(runAllTests ${GROMACS_LIBRARIES} ${GTEST_BOTH_LIBRARIES} pthread)

    # make ctest aware of test executable:
    # (if the color option is not set some test fail)
    add_test(runAllTests runAllTests --gtest_color=yes)

    # also explicitly create binary:
    # TODO: remove this in release:
    install(TARGETS runAllTests DESTINATION /sansom/s117/scro2967/repos/chap/bin)

endif()

