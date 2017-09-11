#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>

#include <gromacs/utility/programcontext.h>

extern const char g_CHAP_INSTALL_BASE[];


/*!
 *
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

