#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>

#include <gromacs/utility/programcontext.h>

extern const char g_CHAP_INSTALL_BASE[];


/*!
 * \brief Returns base directory of install.
 */
extern "C" inline std::string
chapInstallBase()
{
    return std::string(g_CHAP_INSTALL_BASE);
}


/*!
 * \brief Returns command line used to start a given execution of program.
 */
extern "C" inline std::string
chapCommandLine()
{
    // get command line from Gromacs library:
    std::string cmdl = gmx::getProgramContext().commandLine();

    // hack the -quiet flag:
    cmdl.resize(cmdl.size() - 7);   

    return cmdl;
}

#endif

