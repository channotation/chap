#ifndef CHAP_COMMAND_LINE_RUNNER_HPP
#define CHAP_COMMAND_LINE_RUNNER_HPP

#include <functional>
#include <memory>

#include <gromacs/commandline/cmdlineoptionsmodule.h>
#include <gromacs/trajectoryanalysis.h>




/*
 *
 */
class ChapTrajAnaCommandLineRunner
{
    public:

        //
        typedef std::function<gmx::TrajectoryAnalysisModulePointer()> ModuleFactoryMethod;

        //
        template <class ModuleType>
        static int runAsMain(int argc, char *argv[])
        {
            return runAsMain(argc, argv, &createModule<ModuleType>);
        }

        //
        static int runAsMain(int argc, char *argv[], ModuleFactoryMethod factory);

        //
        static std::unique_ptr<gmx::ICommandLineOptionsModule>
        createModule(gmx::TrajectoryAnalysisModulePointer module);

    private:

        // disallow instantiation:
        ChapTrajAnaCommandLineRunner(){};

        //
        template <class ModuleType>
        static gmx::TrajectoryAnalysisModulePointer createModule()
        {
            return TrajectoryAnalysisModulePointer(new ModuleType());
        }
 
};

#endif

