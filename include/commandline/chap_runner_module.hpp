#ifndef CHAP_RUNNER_MODULE_HPP
#define CHAP_RUNNER_MODULE_HPP

#include <gromacs/commandline/cmdlineoptionsmodule.h>
#include <gromacs/trajectoryanalysis.h>

#include "commandline/chap_topology_provider.hpp"
#include "commandline/chap_traj_ana_runner_common.hpp"
#include "trajectory-analysis/trajectory-analysis.hpp"

/*
 *
 */
class ChapRunnerModule : public gmx::ICommandLineOptionsModule
{
    public:

        //
        explicit ChapRunnerModule(ChapTrajectoryAnalysisModulePointer module)
            : module_(std::move(module)), common_(&settings_) {};

        // implementation of interface specified by ICommandLineOptionsModule:
        virtual void init(gmx::CommandLineModuleSettings *settings);
        virtual void initOptions(gmx::IOptionsContainer *options,
                                 gmx::ICommandLineOptionsModuleSettings *settings);
        virtual void optionsFinished();
        virtual int run();
   
        ChapTrajectoryAnalysisModulePointer module_;
        ChapTrajAnaRunnerCommon common_;
        gmx::TrajectoryAnalysisSettings settings_;
        gmx::SelectionCollection selections_;
    

        ChapTopologyProvider topologyProvider_;
 

};

#endif

