#ifndef PDB_IO_HPP
#define PDB_IO_HPP

#include <fstream>
#include <string>
#include <vector>

#include "gromacs/utility/real.h"


/*!
 * \brief Data container for one atom line in a PDB file.
 */
class PdbAtom
{
    friend class PdbIo;
    friend class PdbStructure;

    public:

        // constructor:
        PdbAtom(
                const int serial,
                const std::string atomName,
                const std::string altLoc,
                const std::string resName,
                const std::string chain,
                const int resSeq,
                const std::string iCode,
                const real x,
                const real y,
                const real z,
                const real occupancy,
                const real tempFactor,
                std::string element,
                std::string charge);

    private:

        // atom line attributes:
        int serial_;
        std::string atomName_;
        std::string altLoc_;
        std::string resName_;
        std::string chain_;
        int resSeq_;
        std::string iCode_;
        real x_;
        real y_;
        real z_;
        real occupancy_;
        real tempFactor_;
        std::string element_;
        std::string charge_;

};


/*!
 * \brief Container class for a PDB structure.
 */
class PdbStructure
{
    friend class PdbIo;

    public:

        // setter functions:
        void addAtom(
                const int serial,
                const std::string atomName,
                const std::string altLoc,
                const std::string resName,
                const std::string chain,
                const int resSeq,
                const std::string iCode,
                const real x,
                const real y,
                const real z,
                const real occupancy,
                const real tempFactor,
                std::string element,
                std::string charge);

    private:

        std::vector<PdbAtom> atoms_;
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

        // file handle:
        std::fstream file_;

        // helper functions for writing specific PDB entries:
        inline void writeAtom(const PdbAtom &atom);
        inline void writeRemark(std::string remark);
};

#endif

