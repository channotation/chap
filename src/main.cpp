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


#include <iostream>
#include <string>
#include <vector>

//#include <gromacs/trajectoryanalysis.h>
#include <gromacs/commandline/cmdlineinit.h>
#include <gromacs/commandline/cmdlineprogramcontext.h>
#include <gromacs/commandline.h>

#include "config/version.hpp"
#include "config/config.hpp"
#include "trajectory-analysis/trajectory-analysis.hpp"

using namespace gmx;


int main(int argc, char **argv)
{
    // hack to suppress Gromacs output:
    std::vector<char*> modArgv(argv, argv + argc);
    char quiet[7] = "-quiet";
    modArgv.push_back(quiet);
    modArgv.push_back(nullptr);
    argv = modArgv.data();
    argc++;

    // run trajectory analysis:
	int status =  gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<trajectoryAnalysis>(argc, argv);

    // return status:
    return status;
}











