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


#include <vector>

#include "config/back_matter.hpp"
#include "config/front_matter.hpp"
#include "trajectory-analysis/chap_trajectory_analysis.hpp"

using namespace gmx;


int main(int argc, char **argv)
{
    // print front matter:
    FrontMatter::print();

    // hack to suppress Gromacs output:
    std::vector<char*> modArgv(argv, argv + argc);
    char quiet[7] = "-quiet";
    modArgv.push_back(quiet);
    modArgv.push_back(nullptr);
    argv = modArgv.data();
    argc++;

    // run trajectory analysis:
	int status = TrajectoryAnalysisCommandLineRunner::runAsMain<ChapTrajectoryAnalysis>(argc, argv);

    // print back matter:
    BackMatter::print();

    // return status:
    return status;
}

