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


#include <cstring>
#include <iomanip>

#include <gromacs/fileio/confio.h>
#include <gromacs/topology/atoms.h>

#include "io/pdb_io.hpp"


/*!
 * Creates a PdbStructure from the coordinates contained in the given 
 * TopologyInformation.
 */
void
PdbStructure::fromTopology(
        const gmx::TopologyInformation &top)
{
    // retrieve coordinates and box matrix from topology:
    top.getTopologyConf(&coords_, box_);

    // retrieve list of atoms in topology:
    t_topology *topol = top.topology();
    atoms_ = topol -> atoms;  

    // retrieve periodic BC:
    ePBC_ = top.ePBC();
}


/*!
 * Sets the occupancy and bfac fields of the PDB file to the time averaged 
 * pore-lining and pore-facing attributes.
 */
void
PdbStructure::setPoreFacing(
        const std::vector<SummaryStatistics> &poreLining,
        const std::vector<SummaryStatistics> &poreFacing)
{
    // if no pdbinfo exists, need to create array manually:
    if( atoms_.pdbinfo == nullptr )
    {    
        // allocate array:
        atoms_.pdbinfo = new t_pdbinfo[atoms_.nr];

        // fill with proper values (occup and bfac assigned zero by default):
        for(size_t i = 0; i < atoms_.nr; i++)
        {
            atoms_.pdbinfo[i].type = epdbATOM;
            atoms_.pdbinfo[i].atomnr = i + 1;
            atoms_.pdbinfo[i].altloc = ' ';
            atoms_.pdbinfo[i].occup = 0.0;
            atoms_.pdbinfo[i].bfac = 0.0;
            std::strncpy(atoms_.pdbinfo[i].atomnm, *(atoms_.atomname[i]), 5);
            atoms_.pdbinfo[i].bAnisotropic = false;

        }
    }

    // assign pore facing/lining to occupancy and bfac:
    for(size_t i = 0; i < atoms_.nr; i++)
    {
        // find residue of this atom:
        size_t resind = atoms_.atom[i].resind;

        // have we measured pore facing attribute for this residue?
        if( resind >= poreLining.size() )
        {
            // both attributes zero by default:
            atoms_.pdbinfo[i].occup = 0.0;
            atoms_.pdbinfo[i].bfac = 0.0;
        }
        else
        {
            // reassign occupancy and bfactor:
            atoms_.pdbinfo[i].occup = poreLining.at(resind).mean();
            atoms_.pdbinfo[i].bfac = poreFacing.at(resind).mean();
        }
    }
}


/*!
 * Writes a given PDB structure to a file.
 */
void
PdbIo::write(
        std::string fileName,
        PdbStructure structure)
{
    // delete file if it already exists:
    std::remove(fileName.c_str());

    // create a title string for the PDB file:
    std::string title("created by CHAP");

    // use libgromacs to actually write the PDB file:
    write_sto_conf(
            fileName.c_str(),       // out put file name
            title.c_str(),          // title entry in PDB
            &structure.atoms_,      // atoms in topology
            structure.coords_,      // atom coordinates
            NULL,                   // velocities
            structure.ePBC_,        // periodic BC
            structure.box_);        // box matrix
}

