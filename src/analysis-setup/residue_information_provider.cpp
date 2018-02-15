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


#include "analysis-setup/residue_information_provider.hpp"

#include <cmath>
#include <iostream>


/*!
 * Constructor sets default hydrophobicity to NaN.
 */
ResidueInformationProvider::ResidueInformationProvider()
    : defaultHydrophobicity_(std::nan(""))
{

}


/*!
 * Extracts the residue name for each residue from the topology and adds it to
 * an internal storage container.
 */
void
ResidueInformationProvider::nameFromTopology(
        const gmx::TopologyInformation &top)
{
    // get list of all atoms (including residue information):
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
        "hydrophobicity": [
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
 * Set default hydrophobicity that will be used for residues for which no 
 * lookup table entry exists.
 */
void
ResidueInformationProvider::setDefaultHydrophobicity(const real hydrophobicity)
{
    defaultHydrophobicity_ = hydrophobicity;
}


/*!
 * Returns vector of all IDs for which a name is known.
 */
std::vector<int>
ResidueInformationProvider::ids() const
{
    // loop over name field and add all keys to vector:
    std::vector<int> ids;
    ids.reserve(name_.size());
    for(auto &pair : name_)
    {
        ids.push_back(pair.first);
    }

    // return vector of IDs:
    return(ids);
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
 * Returns hydrophobicity of residue of given ID. This first checks if the 
 * internal lookup table generated with hydrophobicityFromJson() contains an
 * entry for the specified residue. If corresponding value is found, it checks
 * if a fallback hydrophobicity has been specified with 
 * setDefaultHydrophobicity(). If this is not the case, it throws an exception.
 */
real
ResidueInformationProvider::hydrophobicity(const int id) const
{
    // check if record is present:
    if( hydrophobicity_.find(name(id)) == hydrophobicity_.end() )
    {
        // check if default has been set:
        if( !std::isnan(defaultHydrophobicity_) )
        {
            return defaultHydrophobicity_;
        }
        else
        {
            throw std::runtime_error("No hydrophobicity scale data found for "
            "residue " + name(id) + " and no fallback specified.");
        }
    }
    else
    {
        // return lookup value if found:
        return hydrophobicity_.at( name(id) );
    }
}

