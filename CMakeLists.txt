# minimum version of cmake required for building:
cmake_minimum_required(VERSION 2.8.8)

# set project name:
project(CHAP)


# set path to custom cmake modules:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# bring header files into build environment:
include_directories(${CHAP_SOURCE_DIR}/include)

# set library directories for the linker:
link_directories( /usr/local/fftw/lib/ )
#link_directories( /usr/include/ )


# add source files:
#file(GLOB SOURCES "src/*.cpp")

# set compiler and linker flags:
set(CMAKE_CXX_FLAGS "-Wall")

set(LAPACK_LINKER_FLAGS "-llapacke -llapack -lblas")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}" )


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


# Find Linear Algebra Libraries
#------------------------------------------------------------------------------

# find LAPACK:
find_package(LAPACK REQUIRED)

message("lapack = " ${LAPACK_LIBRARIES})


# Compile Instructions
#------------------------------------------------------------------------------

# create executable chap from main.cpp:
add_executable(chap src/main.cpp src/trajectoryAnalysis/trajectoryAnalysis.cpp src/trajectoryAnalysis/simulatedAnnealingModule.cpp)
message("COMPILE_FLAGS = " ${GROMACS_CXX_FLAGS})
set_target_properties(chap PROPERTIES  COMPILE_FLAGS "${GROMACS_CXX_FLAGS} -llapack")
message(STATUS "libs = " ${GROMACS_LIBRARIES})
target_link_libraries(chap ${GROMACS_LIBRARIES})


# Compile Tests
#------------------------------------------------------------------------------

# add test subdirectory:
add_subdirectory(test)


# Install Destinations
#------------------------------------------------------------------------------

# where to install executable on the system:
install(TARGETS chap DESTINATION /sansom/s117/scro2967/repos/chap/bin)

