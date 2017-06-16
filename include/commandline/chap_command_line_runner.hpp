#ifndef CHAP_COMMAND_LINE_RUNNER_HPP
#define CHAP_COMMAND_LINE_RUNNER_HPP

#include <functional>
#include <memory>

#include <gromacs/commandline/cmdlineoptionsmodule.h>
#include <gromacs/trajectoryanalysis.h>

#include "trajectory-analysis/trajectory-analysis.hpp"



/*
 *
 */
class ChapTrajAnaCommandLineRunner
{
    public:

        //
        typedef std::function<ChapTrajectoryAnalysisModulePointer()> ModuleFactoryMethod;

        //
        template <class ModuleType>
        static int runAsMain(int argc, char *argv[])
        {
            return runAsMain(argc, argv, &createModule<ModuleType>);
        }

        //
        static int runAsMain(int argc, char *argv[], ModuleFactoryMethod factory);

        //
        static void registerModule(gmx::CommandLineModuleManager *manager,
                                   const char *name,
                                   const char *description,
                                   ModuleFactoryMethod factory);

        //
        static std::unique_ptr<gmx::ICommandLineOptionsModule>
        createModule(ChapTrajectoryAnalysisModulePointer module);

    private:

        // disallow instantiation:
        ChapTrajAnaCommandLineRunner(){};

        //
        template <class ModuleType>
        static ChapTrajectoryAnalysisModulePointer createModule()
        {
            return ChapTrajectoryAnalysisModulePointer(new ModuleType());
        }
 
};

#endif

