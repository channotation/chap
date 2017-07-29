# minimum version of cmake required for building:
cmake_minimum_required(VERSION 2.8.8)

# set project name:
project(CHAP)


# Prepare Version Header and Config Header
#------------------------------------------------------------------------------

# set path to custom cmake modules:
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

# include cmake scripts for geting git revisions:
include(GetGitRevisionDescription)

# version definitions:
set(CHAP_VERSION_MAJOR "0")
set(CHAP_VERSION_MINOR "2")
set(CHAP_VERSION_PATCH "0")

# get git hash:
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# configure header files:
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/config/version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/config/version.cpp" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/config/config.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/config/config.cpp" @ONLY)


# Find Linear Algebra Libraries
#------------------------------------------------------------------------------

# find LAPACK:
find_package(LAPACKE REQUIRED)
include_directories(${LAPACKE_INCLUDE_DIRS})


# Find Boost Library
#------------------------------------------------------------------------------

find_package(Boost REQUIRED)


# Find Google Test Library
#------------------------------------------------------------------------------

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})


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
find_package(GROMACS 5.1 REQUIRED)
gromacs_check_double(GMX_DOUBLE)
gromacs_check_compiler(CXX)
include_directories(${GROMACS_INCLUDE_DIRS})
add_definitions(${GROMACS_DEFINITIONS})

# 
#------------------------------------------------------------------------------


# bring header files into build environment:
include_directories(${CHAP_SOURCE_DIR}/include)

# set library directories for the linker:
link_directories( /usr/local/fftw/lib/ )
#link_directories( /usr/include/ )


# add source files:
#file(GLOB SOURCES "src/*.cpp")

# set compiler and linker flags:
set(CMAKE_CXX_FLAGS "-Wall -Wextra -O1 -g")

set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}" )





# Compile Instructions
#------------------------------------------------------------------------------

# rapidjson support for std::string:
add_definitions(-DRAPIDJSON_HAS_STDSTRING)

file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
list(APPEND SRC_FILES "${CMAKE_CURRENT_BINARY_DIR}/config/version.cpp")
list(APPEND SRC_FILES "${CMAKE_CURRENT_BINARY_DIR}/config/config.cpp")
#message("CHAP_SRC = " ${SRC_FILES})

# create executable chap from main.cpp:
add_executable(chap ${SRC_FILES})
#message("COMPILE_FLAGS = " ${GROMACS_CXX_FLAGS})
set_target_properties(chap PROPERTIES  COMPILE_FLAGS "${GROMACS_CXX_FLAGS} -llapack")
#message(STATUS "libs = " ${GROMACS_LIBRARIES})
target_link_libraries(chap ${LAPACKE_LIBRARIES})
target_link_libraries(chap ${BOOST_LIBRARIES})
target_link_libraries(chap ${GTEST_BOTH_LIBRARIES})
target_link_libraries(chap ${GROMACS_LIBRARIES})


# Compile Tests
#------------------------------------------------------------------------------

# include(CTest) rather than enable_testing() to allow calls to make test:
# (this does not work if called from CMakeLists.txt in test subdirectory)
include(CTest)

# add test subdirectory:S
add_subdirectory(test)


# Install Destinations
#------------------------------------------------------------------------------

# TODO: remove this for release
#set(CMAKE_INSTALL_PREFIX "/sansom/s117/scro2967/repos/chap/bin")

# where to install executable on the system:
install(TARGETS chap DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# also install data:
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/share DESTINATION ${CMAKE_INSTALL_PREFIX})

