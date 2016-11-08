#include <iostream>
#include <string>
#include <vector>

#include <gromacs/trajectoryanalysis.h>

#include "trajectoryAnalysis/trajectoryAnalysis.hpp"

using namespace gmx;








int main(int argc, char **argv)
{
	std::cout<<"Hello, beautiful world!"<<std::endl;

	int status =  gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<trajectoryAnalysis>(argc, argv);

	std::cout<<"done. status = "<<status<<std::endl;

	return 0;
}











