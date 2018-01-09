#include <iomanip>

#include <gromacs/fileio/confio.h>
#include <gromacs/topology/atoms.h>

#include "io/pdb_io.hpp"


// FIXME remove:
#include <iostream>


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
 * Writes a given PDB structure to a file.
 */
void
PdbIo::write(
        std::string fileName,
        PdbStructure structure)
{
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

