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
#include "commandline/chap_runner_module.hpp"
#include "commandline/chap_command_line_runner.hpp"



using namespace gmx;





int main(int argc, char **argv)
{
	std::cout<<"Hello, beautiful world!"<<std::endl;


//    gmx::CommandLineProgramContext &programContext = gmx::initForCommandLine(&argc, &argv);


//    ChapRunnerModule();


    

    

//	int status =  gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<trajectoryAnalysis>(argc, argv);


	int status =  ChapTrajAnaCommandLineRunner::runAsMain<trajectoryAnalysis>(argc, argv);


/*
    gmx::CommandLineProgramContext &programContext = gmx::initForCommandLine(&argc, &argv);

    try
    {
        gmx::CommandLineModuleManager manager(NULL, &programContext);
//        int rc = manager.run();
        gmx::finalizeForCommandLine();
//        return rc;
    }
    catch(const std::exception &ex)
    {
        gmx::printFatalErrorMessage(stderr, ex);
        return gmx::processExceptionAtExitForCommandLine(ex);
    }
  
*/


/*

    std::cout<<"SHA1 = "<<g_GIT_SHA1<<std::endl;
    std::cout<<"CHAP_INSTALL_BASE = "<<g_CHAP_INSTALL_BASE<<std::endl;


    CommandLineProgramContext &test = gmx::initForCommandLine(&argc, &argv);

    std::cout<<"full binary path = "<<std::string(test.fullBinaryPath())<<std::endl;
    std::cout<<"program name = "<<std::string(test.programName())<<std::endl;
    std::cout<<"display name = "<<std::string(test.displayName())<<std::endl;

    gmx::finalizeForCommandLine();

	int status =  gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<trajectoryAnalysis>(argc, argv);

	std::cout<<"done. status = "<<status<<std::endl;

	return 0;
    */
}











