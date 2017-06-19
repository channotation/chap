#include <iostream>
#include <string>
#include <vector>

//#include <gromacs/trajectoryanalysis.h>
#include <gromacs/commandline/cmdlineinit.h>
#include <gromacs/commandline/cmdlineprogramcontext.h>
#include <gromacs/commandline.h>

#include "config/version.hpp"
#include "config/config.hpp"
#include "trajectory-analysis/trajectory-analysis.hpp"

using namespace gmx;


int main(int argc, char **argv)
{
	int status =  gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<trajectoryAnalysis>(argc, argv);
}











