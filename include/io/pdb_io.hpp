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


#ifndef PDB_IO_HPP
#define PDB_IO_HPP

#include <fstream>
#include <string>
#include <vector>

#include <gromacs/trajectoryanalysis/analysissettings.h>
#include <gromacs/topology/topology.h>
#include <gromacs/utility/real.h>

#include "statistics/summary_statistics.hpp"


/*!
 * \brief Container class for a PDB structure.
 */
class PdbStructure
{
    friend class PdbIo;

    public:

        // create PDB file from topology:
        void fromTopology(const gmx::TopologyInformation &top);

        // 
        void setPoreFacing(
                const std::vector<SummaryStatistics> &poreLining,
                const std::vector<SummaryStatistics> &poreFacing);


    private:

        // data required for writing PDB file:
        t_atoms atoms_;
        rvec *coords_;
        int ePBC_;
        matrix box_;
};


/*!
 * \brief Exports structures to PDB file format.
 *
 * This wraps around the PDB export utilities of Gromacs.
 */
class PdbIo
{
    public:

        // public interface for PDB export:
        static void write(
                std::string fileName,
                PdbStructure structure);

    private:

};

#endif

