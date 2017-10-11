#include "analysis-setup/residue_information_provider.hpp"


#include <iostream>


/*!
 * Extracts the residue name for each residue from the topology and adds it to
 * an internal storage container.
 */
void
ResidueInformationProvider::nameFromTopology(
        const gmx::TopologyInformation &top)
{
    // get list of all atoms (includingb residue information):
    t_atoms atoms = top.topology() -> atoms;

    // loop over all residues:
    for(int i = 0; i < 800; i++)
    {
        // add chain ID to map:
        name_[i] = std::string(*atoms.resinfo[i].name);
    }
}


/*!
 * Extracts chain ID for each residue from the topology and adds it to an
 * internal storage container.
 */
void
ResidueInformationProvider::chainFromTopology(
        const gmx::TopologyInformation &top)
{
    // get list of all atoms (includingb residue information):
    t_atoms atoms = top.topology() -> atoms;

    // loop over all residues:
    for(int i = 0; i < atoms.nres; i++)
    {
        // add chain ID to map:
        chain_[i] = atoms.resinfo[i].chainid;
    }
}


/*
 *
 */
void
ResidueInformationProvider::hydrophobicityFromJson(
        const rapidjson::Document &doc)
{
    
}


/*!
 * Returns name of residue of given ID.
 */
std::string
ResidueInformationProvider::name(const int id) const
{
    return name_.at(id);
}


/*!
 * Returns chain ID of residue of given ID.
 */
std::string
ResidueInformationProvider::chain(const int id) const
{
    return chain_.at(id);
}


/*!
 * Returns hydrophobicity of residue of gievn ID.
 */
real
ResidueInformationProvider::hydrophobicity(const int id) const
{
    return hydrophobicity_.at(id);
}

