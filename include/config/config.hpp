// CHAP - The Channel Annotation Package
// 
// Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
// Stephen J. Tucker
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


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

