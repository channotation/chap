# CHAP - The Channel Annotation Package
# 
# Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
# Stephen J. Tucker
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


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

# need pthreads for Google test:
find_package(Threads)

# add executable to run all tests and link libraries:
add_executable(runAllTests ${TEST_SRC_FILES} ${SRC_FILES})
target_include_directories(runAllTests PUBLIC ${CHAP_SOURCE_DIR}/include)
target_link_libraries(runAllTests ${GROMACS_LIBRARIES})
target_link_libraries(runAllTests ${LAPACKE_LIBRARIES})
target_link_libraries(runAllTests ${LAPACK_LIBRARIES})
target_link_libraries(runAllTests ${BLAS_LIBRARIES})
target_link_libraries(runAllTests ${GTEST_LIBRARY})
target_link_libraries(runAllTests ${CMAKE_THREAD_LIBS_INIT})

# make ctest aware of test executable:
# (if the color option is not set some test fail)
add_test(runAllTests runAllTests --gtest_color=yes)

