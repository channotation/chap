#include "commandline/chap_command_line_runner.hpp"
#include "commandline/chap_runner_module.hpp"


// TODO: this is suicide
//#include "gromacs/"


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
void
ChapTrajAnaCommandLineRunner::registerModule(
        gmx::CommandLineModuleManager *manager,
        const char *name,
        const char *description,
        ModuleFactoryMethod factory)
{
    auto runnerFactory = [factory]
    {
        return createModule(factory());
    };
    ICommandLineOptionsModule::registerModuleFactory(manager,
                                                     name,
                                                     description,
                                                     runnerFactory);
}


/*
 *
 */
std::unique_ptr<gmx::ICommandLineOptionsModule>
ChapTrajAnaCommandLineRunner::createModule(ChapTrajectoryAnalysisModulePointer module)
{
    return gmx::ICommandLineOptionsModulePointer(new ChapRunnerModule(std::move(module)));
}

