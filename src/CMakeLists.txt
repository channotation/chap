# minimum version of cmake required for building:
cmake_minimum_required(VERSION 2.6)

# set project name:
project(CHAP)

# set compiler flags:
set(CMAKE_CXX_FLAGS "-Wall")

# create executable chap from main.cpp:
add_executable(chap main.cpp)

# where to install executable on the system:
install(TARGETS chap DESTINATION /sansom/s117/scro2967/repos/chap/bin)
