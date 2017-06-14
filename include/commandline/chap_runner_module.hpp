#ifndef CHAP_RUNNER_MODULE_HPP
#define CHAP_RUNNER_MODULE_HPP

#include <gromacs/commandline/cmdlineoptionsmodule.h>
#include <gromacs/trajectoryanalysis.h>



/*
 *
 */
class ChapRunnerModule : public gmx::ICommandLineOptionsModule
{
    public:

        //
        explicit ChapRunnerModule(gmx::TrajectoryAnalysisModulePointer module)
            : module_(std::move(module)) {};

        // implementation of interface specified by ICommandLineOptionsModule:
        virtual void init(gmx::CommandLineModuleSettings *settings);
        virtual void initOptions(gmx::IOptionsContainer *options,
                                 gmx::ICommandLineOptionsModuleSettings *settings);
        virtual void optionsFinished();
        virtual int run();
   
        gmx::TrajectoryAnalysisModulePointer module_;
        gmx::TrajectoryAnalysisSettings settings_;
//        extern TrajectoryAnalysisRunnerCommon common_;
        gmx::SelectionCollection selections_;
    

};

#endif

