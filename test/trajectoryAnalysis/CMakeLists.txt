# unit tests for the simulated annealing module:
add_executable(utSimulatedAnnealingModule utSimulatedAnnealingModule.cpp)
target_link_libraries(utSimulatedAnnealingModule ${GTEST_LIBRARIES} pthread)
install(TARGETS utSimulatedAnnealingModule DESTINATION /sansom/s117/scro2967/repos/chap/bin)
