# unit tests for the simulated annealing module:
add_executable(utSimulatedAnnealingModule utSimulatedAnnealingModule.cpp ${CHAP_SOURCE_DIR}/src/trajectoryAnalysis/simulatedAnnealingModule.cpp)
target_link_libraries(utSimulatedAnnealingModule ${GROMACS_LIBRARIES} ${GTEST_LIBRARIES} pthread)
install(TARGETS utSimulatedAnnealingModule DESTINATION /sansom/s117/scro2967/repos/chap/bin)
