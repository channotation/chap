#include <iostream>
#include <string>
#include <vector>

#include <gromacs/trajectoryanalysis.h>

#include "config/version.hpp"
#include "trajectory-analysis/trajectory-analysis.hpp"

using namespace gmx;








int main(int argc, char **argv)
{
	std::cout<<"Hello, beautiful world!"<<std::endl;


    std::cout<<"SHA1 = "<<g_GIT_SHA1<<std::endl;


	int status =  gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<trajectoryAnalysis>(argc, argv);

	std::cout<<"done. status = "<<status<<std::endl;

	return 0;
}











