# Setup Google Test
# -----------------------------------------------------------------------------

# create target for running make check:
add_custom_target(check ${CMAKE_CTEST_COMMAND} -V)

# get list of all source files but ignore main file of chap directory:
file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB_RECURSE TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
list(REMOVE_ITEM SRC_FILES ${PROJECT_SOURCE_DIR}/src/main.cpp)
list(APPEND SRC_FILES "${CMAKE_CURRENT_BINARY_DIR}/../config/version.cpp")
list(APPEND SRC_FILES "${CMAKE_CURRENT_BINARY_DIR}/../config/config.cpp")

# add executable to run all tests and link libraries:
add_executable(runAllTests ${TEST_SRC_FILES} ${SRC_FILES})
set_target_properties(runAllTests PROPERTIES COMPILE_FLAGS ${GROMACS_CXX_FLAGS})
target_include_directories(runAllTests PUBLIC ${CHAP_SOURCE_DIR}/include)
target_link_libraries(runAllTests ${GROMACS_LIBRARIES} pthread)
target_link_libraries(runAllTests ${LAPACKE_LIBRARIES})
target_link_libraries(runAllTests ${LAPACK_LIBRARIES})
target_link_libraries(runAllTests ${BLAS_LIBRARIES})
target_link_libraries(runAllTests ${GTEST_LIBRARY})

# make ctest aware of test executable:
# (if the color option is not set some test fail)
add_test(runAllTests runAllTests --gtest_color=yes)

