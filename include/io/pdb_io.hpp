#ifndef PDB_IO_HPP
#define PDB_IO_HPP

#include <fstream>
#include <string>
#include <vector>

#include <gromacs/trajectoryanalysis/analysissettings.h>
#include <gromacs/trajectory/trajectoryframe.h>
#include <gromacs/topology/topology.h>
#include <gromacs/utility/real.h>



/*!
 * \brief Container class for a PDB structure.
 */
class PdbStructure
{
    friend class PdbIo;

    public:

        // create PDB file from topology:
        void fromTopology(const gmx::TopologyInformation &top);


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
 * This functionality is of course also part of Gromacs, but not of the public
 * library API. To avoid that Gromacs changes the interface, this is my own 
 * simplistic implementation. 
 */
class PdbIo
{
    public:

        // public interface for PDB export:
        void write(
                std::string fileName,
                PdbStructure structure);

    private:

};

#endif

