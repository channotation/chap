#include <vector>

#include "config/back_matter.hpp"
#include "config/front_matter.hpp"
#include "trajectory-analysis/chap_trajectory_analysis.hpp"

using namespace gmx;


int main(int argc, char **argv)
{
    // print front matter:
    FrontMatter::print();

    // hack to suppress Gromacs output:
    std::vector<char*> modArgv(argv, argv + argc);
    char quiet[7] = "-quiet";
    modArgv.push_back(quiet);
    modArgv.push_back(nullptr);
    argv = modArgv.data();
    argc++;

    // run trajectory analysis:
	int status = TrajectoryAnalysisCommandLineRunner::runAsMain<ChapTrajectoryAnalysis>(argc, argv);

    // print back matter:
    BackMatter::print();

    // return status:
    return status;
}

