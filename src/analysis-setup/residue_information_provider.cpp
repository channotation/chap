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
    for(int i = 0; i < atoms.nres; i++)
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


/*!
 * Creates an internal hydrophobicity lookup table from a JSON document. This 
 * table needs to provide a hydrophobicity value for each residue in the pore
 * forming molecule (typically one for each amino acid) and may not contain 
 * duplicate entries. The layout of the JSON document required as input is as
 * follows:
 *
 *  \code{.java} 
    {
        "vdwradii": [
            {"resname": "ALA", "hydrophobicity":  41.0},
            {"resname": "ARG", "hydrophobicity": -14.0}
        ]
    }
 *  \endcode
 *
 * The order of the records is irrelevant.
 */
void
ResidueInformationProvider::hydrophobicityFromJson(
        const rapidjson::Document &doc)
{
    // sanity checks:
    if( !doc.IsObject() )
    {
        throw std::runtime_error("No valid JSON object provide for generation "
        "of hydrophobicity scale.");
    }
    if( !doc.HasMember("hydrophobicity") || !doc["hydrophobicity"].IsArray() )
    {
        throw std::runtime_error("JSON document provided for hydrophobicity "
        "scale generation does not contain hydrophobicity array.");
    }

    // extract hydrophobicity data:
    const rapidjson::Value &hydrophobicityEntries = doc["hydrophobicity"];

    // iterate over provided values:
    rapidjson::Value::ConstValueIterator it;
    for(it = hydrophobicityEntries.Begin(); 
        it != hydrophobicityEntries.End();
        it++)
    {
        // check that required entries are present and of correct type:
        if( !(it -> HasMember("resname")) || 
            !(*it)["resname"].IsString() )
        {
            throw std::runtime_error("No 'resname' attribute of type string "
            "in hydrophobicity record");
        }
        if( !(it -> HasMember("hydrophobicity")) || 
            !(*it)["hydrophobicity"].IsDouble() )
        {
            throw std::runtime_error("No 'resname' attribute of type string "
            "in hydrophobicity record");
        }

        // prevent duplicate entries:
        std::string resname = (*it)["resname"].GetString();
        if( hydrophobicity_.find(resname) == hydrophobicity_.end() )
        {
            // add to internal lookup table:
            hydrophobicity_[resname] = (*it)["hydrophobicity"].GetDouble(); 
        }
        else
        {
            throw std::runtime_error("Duplicate entry in hydrophobicity "
            "database.");
        }
    }
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
    return hydrophobicity_.at( name(id) );
}

