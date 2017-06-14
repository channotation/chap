#include "commandline/chap_command_line_runner.hpp"
#include "commandline/chap_runner_module.hpp"

// TODO: this is suicide
//#include "gromacs/"

/*
 *
 */
/*
template <class ModuleType>
int
ChapTrajAnaCommandLineRunner::runAsMain(int argc, 
                                        char *argv[])
{
    return runAsMain(argc, argv, &createModule<ModuleType>);
}
*/

/*
 *
 */
int
ChapTrajAnaCommandLineRunner::runAsMain(int argc, 
                                        char *argv[],
                                        ModuleFactoryMethod factory)
{
    auto runnerFactory = [factory]
    {
        return createModule(factory());
    };

    return gmx::ICommandLineOptionsModule::runAsMain(argc, argv, NULL, NULL, runnerFactory);
}


/*
 *
 */
std::unique_ptr<gmx::ICommandLineOptionsModule>
ChapTrajAnaCommandLineRunner::createModule(gmx::TrajectoryAnalysisModulePointer module)
{
    return gmx::ICommandLineOptionsModulePointer(new ChapRunnerModule(std::move(module)));
}


/*
 *
 */
/*
template <class ModuleType>
gmx::TrajectoryAnalysisModulePointer
ChapTrajAnaCommandLineRunner::createModule()
{
    return TrajectoryAnalysisModulePointer(new ModuleType());
}
*/

