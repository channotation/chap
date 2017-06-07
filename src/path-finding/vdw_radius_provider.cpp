#include <iostream>

//#include <gromacs/topology/atomprop.h> 
#include <gromacs/topology/atoms.h>  
#include <gromacs/topology/topology.h>  

#include "path-finding/vdw_radius_provider.hpp"


/*
 *
 */
VdwRadiusProvider::VdwRadiusProvider()
    : defRad_(-1.0)
{

}

/*
 *
 */
VdwRadiusProvider::~VdwRadiusProvider()
{

}


/*
 *
 */
void
VdwRadiusProvider::lookupTableFromJson(rapidjson::Document &jsonDoc)
{
    // ensure that root of JSON is object:
    if( jsonDoc.IsObject() == false )
    {
        throw std::runtime_error("No valid JSON object provided.");
    }

    // ensure that JSON contains vdwradii array:
    if( jsonDoc.HasMember("vdwradii") == false || jsonDoc["vdwradii"].IsArray() == false )
    {
        throw std::runtime_error("Provided JSON does not contain vdwradii array.");
    }

    // extract vdw radius data:
    const rapidjson::Value &vdwRadiiEntries = jsonDoc["vdwradii"];

    // prepare loopup table for vdW radii:
    vdwRadiusLookupTable_.clear();
//    vdwRadiusLookupTable_.reserve(vdwRadiusEntries.Size());

    // loop over array and extract vdW radius lookup table:
    rapidjson::Value::ConstValueIterator it;
    for(it = vdwRadiiEntries.Begin(); it != vdwRadiiEntries.End(); it++)
    {
        // TODO make this exceptions

        // check that required entries are present and have correct type:
        if( it -> HasMember("atomname") == false )
        {
            std::cerr<<"ERROR: Van der Waals radius record invalid."<<std::endl;
            std::cerr<<"Attribute 'atomname' not found."<<std::endl;
            std::abort();
        }
        if( (*it)["atomname"].IsString() == false )
        {
            std::cerr<<"ERROR: Van der Waals radius record invalid."<<std::endl;
            std::cerr<<"Attribute 'atomname' must be of type string."<<std::endl;
            std::abort();
        }
        if( it -> HasMember("resname") == false )
        {
            std::cerr<<"ERROR: Van der Waals radius record invalid."<<std::endl;
            std::cerr<<"Attribute 'resname' not found."<<std::endl;
            std::abort();
        }
        if( (*it)["resname"].IsString() == false )
        {
            std::cerr<<"ERROR: van der Waals radius record invalid."<<std::endl;
            std::cerr<<"Attribute 'resname' must be of type string.."<<std::endl;
            std::abort();
        }
        if( it -> HasMember("vdwr") == false )
        {
            std::cerr<<"ERROR: Van der Waals radius record invalid."<<std::endl;
            std::cerr<<"Attribute 'vdwr' not found."<<std::endl;
            std::abort();
        }
        if( (*it)["vdwr"].IsNumber() == false )
        {
            std::cerr<<"ERROR: Van der Waals radius record invalid."<<std::endl;
            std::cerr<<"Attribute 'vdwr' must be of type number."<<std::endl;
            std::abort();
        }

        // create a vdW radius record struct from JSON:
        VdwRadiusRecord rec;
        rec.atmName_ = (*it)["atomname"].GetString();
        rec.resName_ = (*it)["resname"].GetString();
        rec.vdwRad_ = (*it)["vdwr"].GetDouble();

        // add record to lookup table:
        vdwRadiusLookupTable_.push_back(rec);
    } 
}


/*
 *
 */
std::unordered_map<int, real>
VdwRadiusProvider::vdwRadiiForTopology(const gmx::TopologyInformation &top)
{
    // get list of all atoms:
    t_atoms atoms = top.topology() -> atoms;

    // allocate memory for results vector:
    std::unordered_map<int, real> vdwRadii;
    vdwRadii.reserve(atoms.nr);

    // loop over all atoms in topology and find vdW radii:
    for(int i = 0; i < atoms.nr; i++)
    {
        // get atom and residue names:
        std::string atmName = *(atoms.atomname[i]);
        std::string resName = *(atoms.resinfo[atoms.atom[i].resind].name);

        // add radius for this atom to results vector:
        vdwRadii[i] = vdwRadiusForAtom(atmName, resName);
    }

    // return vector of vdW radii:
    return vdwRadii;
}


/*
 *
 */
real
VdwRadiusProvider::vdwRadiusForAtom(std::string atmName, std::string resName)
{
    // build vector of atom name matches:
    std::vector<VdwRadiusRecord> matches;
    std::vector<VdwRadiusRecord>::iterator it;
    for(it = vdwRadiusLookupTable_.begin(); it != vdwRadiusLookupTable_.end(); it++)
    {
        if(it -> atmName_ == atmName )
        {
            matches.push_back(*it);
        }
    }
   
    // have any atom name matches been found:
    if( matches.size() > 0 )
    {
        // find rasidue name match:
        std::vector<VdwRadiusRecord>::iterator jt;
        for(jt = matches.begin(); jt != matches.end(); jt++)
        {
            if( jt -> resName_ == resName )
            {
                return( jt -> vdwRad_);
            }
        }

        // handle case where no match was found:
        if( jt == matches.end() )
        {
            return returnDefaultRadius(atmName, resName);
        }
    }
    else
    {
        return returnDefaultRadius(atmName, resName);
    }
}


/*
 *
 */
real
VdwRadiusProvider::returnDefaultRadius(std::string atmName, std::string resName)
{
    // check if default radius has been set:
    if( defRad_ >= 0 )
    {
        // just return the default radius:
        return(defRad_);
    }
    else
    {
        std::cerr<<"ERROR: Could not find van der Waals radius for atom "
                 <<"with atom name "<<atmName<<" and residue name "
                 <<resName<<" and default radius is not set."<<std::endl;
        std::abort();
    } 
}

















