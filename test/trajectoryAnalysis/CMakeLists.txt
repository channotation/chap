# unit tests for the simulated annealing module:
add_executable(runTests utSimulatedAnnealingModule.cpp)
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)
message(STATUS "utSimulatedAnnealingModule")
install(TARGETS runTests DESTINATION /sansom/s117/scro2967/repos/chap/bin)
