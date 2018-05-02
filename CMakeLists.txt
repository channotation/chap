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


# minimum version of cmake required for building:
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# set project parameters:
project(CHAP LANGUAGES CXX C VERSION 0.9.1)

# follow GNU guideline on install directories:
include(GNUInstallDirs)

# set C++11 standard and turn of GNU extensions:
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# build type defaults to release:
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

# add custom maintainer mode for compiler and linker flags:
SET(CMAKE_CXX_FLAGS_MAINTAINER "-Wall -Wextra -Werror -Wfatal-errors -pedantic")
SET(CMAKE_C_FLAGS_MAINTAINER "-Wall -Wextra -pedantic")
SET(CMAKE_EXE_LINKER_FLAGS_MAINTAINER "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
SET(CMAKE_SHARED_LINKER_FLAGS_MAINTAINER "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}")

# set path to custom cmake modules:
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindGROMACS")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindLAPACKE")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/GetGitRevisionDescription")


# Find Boost Library and Include it as Imported Target
#------------------------------------------------------------------------------

find_package(Boost REQUIRED)
add_library(boost INTERFACE IMPORTED)
set_property(TARGET boost PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})


# Find Linear Algebra Libraries
#------------------------------------------------------------------------------

# find LAPACK:
find_package(LAPACKE REQUIRED)


# Find Gromacs Library
#------------------------------------------------------------------------------

# handle different Gromacs suffixes:
option(GMX_DOUBLE "Use double precision" OFF)
set(GMX_SUFFIX "" CACHE STRING "Suffix for the GROMACS installation to use (empty for default)")
if (GMX_DOUBLE AND NOT GMX_SUFFIX)
	set(GROMACS_SUFFIX "_d")
else()
	set(GROMACS_SUFFIX ${GMX_SUFFIX})
endif()

# find Gromacs package, add include direcotries and definitions:
find_package(GROMACS 2016 REQUIRED)
gromacs_check_double(GMX_DOUBLE)
gromacs_check_compiler(CXX)
include_directories(${GROMACS_INCLUDE_DIRS})
add_definitions(${GROMACS_DEFINITIONS})


# Build Google Test Library as an External Project
#------------------------------------------------------------------------------

# enable external projects:
include(ExternalProject)

# Google test as external project:
ExternalProject_Add(
    googletest
    URL https://github.com/google/googletest/archive/release-1.7.0.zip
    # Disable install step
    INSTALL_COMMAND ""
)

# get source and binary location of Google test libraries:
ExternalProject_Get_Property(googletest source_dir binary_dir)

# set include and library path variables:
set(GTEST_INCLUDE_DIR ${source_dir}/include)
set(GTEST_LIBRARY_PATH ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a)
set(GTEST_LIBRARY gtest)
include_directories(${GTEST_INCLUDE_DIR})

# make library an install target, define properties and dependencies:
add_library(${GTEST_LIBRARY} UNKNOWN IMPORTED)
set_property(TARGET ${GTEST_LIBRARY} PROPERTY IMPORTED_LOCATION
                ${GTEST_LIBRARY_PATH} )
add_dependencies(${GTEST_LIBRARY} googletest)


# Prepare Version Header and Config Header
#------------------------------------------------------------------------------

# include cmake scripts for geting git revisions:
include(GetGitRevisionDescription)

# version definitions:
set(CHAP_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CHAP_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CHAP_VERSION_PATCH ${PROJECT_VERSION_PATCH})

# get git hash:
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# configure header files:
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/config/version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/config/version.cpp" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/config/config.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/config/config.cpp" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/config/dependencies.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/config/dependencies.hpp" @ONLY)
include_directories(${CMAKE_CURRENT_BINARY_DIR})


# Compile Instructions
#------------------------------------------------------------------------------

# rapidjson support for std::string:
add_definitions(-DRAPIDJSON_HAS_STDSTRING)

# build list of sources:
file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
list(APPEND SRC_FILES "${CMAKE_CURRENT_BINARY_DIR}/config/version.cpp")
list(APPEND SRC_FILES "${CMAKE_CURRENT_BINARY_DIR}/config/config.cpp")

# create executable chap from main.cpp:
add_executable(chap ${SRC_FILES})
target_include_directories(chap PUBLIC ${CHAP_SOURCE_DIR}/include)
target_link_libraries(chap ${LAPACKE_LIBRARIES})
target_link_libraries(chap ${BOOST_LIBRARIES})
target_link_libraries(chap ${GROMACS_LIBRARIES})
target_link_libraries(chap ${GTEST_LIBRARY})


# Compile Tests
#------------------------------------------------------------------------------

# include(CTest) rather than enable_testing() to allow calls to make test:
# (this does not work if called from CMakeLists.txt in test subdirectory)
include(CTest)

# add test subdirectory:S
add_subdirectory(test)


# Install Destinations
#------------------------------------------------------------------------------

# where to install executable on the system:
install(TARGETS chap DESTINATION ${CMAKE_INSTALL_PREFIX}/chap/bin)

# also install data and scripts:
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/share DESTINATION ${CMAKE_INSTALL_PREFIX}/chap)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts DESTINATION ${CMAKE_INSTALL_PREFIX}/chap USE_SOURCE_PERMISSIONS)

