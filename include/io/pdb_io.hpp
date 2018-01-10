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

