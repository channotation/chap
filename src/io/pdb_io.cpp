#include <iomanip>

#include "io/pdb_io.hpp"


/*!
 * Create an atom line record from the given attributes.
 */
PdbAtom::PdbAtom(
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
        std::string charge)
    : serial_(serial)
    , atomName_(atomName)
    , altLoc_(altLoc)
    , resName_(resName)
    , chain_(chain)
    , resSeq_(resSeq)
    , iCode_(iCode)
    , x_(x)
    , y_(y)
    , z_(z)
    , occupancy_(occupancy)
    , tempFactor_(tempFactor)
    , element_(element)
    , charge_(charge)
{

}


/*!
 * Adds an atom to the PDB structure.
 */
void
PdbStructure::addAtom(
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
        std::string charge)
{
    // create an atom record container:
    PdbAtom atom(
            serial,
            atomName,
            altLoc,
            resName,
            chain,
            resSeq,
            iCode,
            x,
            y,
            z,
            occupancy,
            tempFactor,
            element,
            charge);

    // add to vector of atom records:
    atoms_.push_back(atom);
}


/*!
 * Writes a given PDB structure to a file.
 */
void
PdbIo::write(
        std::string fileName,
        PdbStructure structure)
{
    // open file stream for writing:
    file_.open(fileName.c_str(), std::fstream::out);

    writeRemark("produced by CHAP");

    // write atom records:
    for(auto atom : structure.atoms_)
    {
        writeAtom(atom);
    }

    // close file stream again:
    file_.close();
}


/*! 
 * Writes an atom line to a PDB file.
 */
void
PdbIo::writeAtom(const PdbAtom &atom)
{
    file_<<std::setw(6)<<"ATOM  "           // record name
         <<std::setw(4)<<atom.serial_       // serial number
         <<std::setw(1)<<""
         <<std::setw(4)<<atom.atomName_     // atom name
         <<std::setw(1)<<atom.altLoc_       // alternate location name
         <<std::setw(3)<<atom.resName_      // residue name
         <<std::setw(1)<<""
         <<std::setw(1)<<atom.chain_        // chain identifier
         <<std::setw(4)<<atom.resSeq_       // residue sequence number
         <<std::setw(1)<<atom.iCode_        // code for insertion of residues
         <<std::setw(3)<<""
         <<std::setw(8)<<atom.x_            // x-coordinate
         <<std::setw(8)<<atom.y_            // y-coordinate
         <<std::setw(8)<<atom.z_            // z-coordinate
         <<std::setw(5)<<atom.occupancy_    // occupancy
         <<std::setw(5)<<atom.tempFactor_   // temperature factor
         <<std::setw(11)<<""
         <<std::setw(2)<<atom.element_      // element symbol, right justified
         <<std::setw(2)<<atom.charge_       // charge on atom
         <<std::endl;
}


/*!
 * Write a remark entry to file.
 */
void
PdbIo::writeRemark(std::string remark)
{
    file_<<"REMARK "
         <<remark
         <<std::endl;
}

