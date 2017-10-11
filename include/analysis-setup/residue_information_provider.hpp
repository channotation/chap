#ifndef RESIDUE_INFORMATION_PROVIDER_HPP
#define RESIDUE_INFORMATION_PROVIDER_HPP

#include <map>
#include <string>

#include <gromacs/topology/atoms.h>
#include <gromacs/topology/topology.h>
#include <gromacs/trajectoryanalysis/analysissettings.h>
#include <gromacs/utility/real.h>

#include "rapidjson/Document"


/*!
 *
 */
class ResidueInformationProvider
{
    public:

        // setter methods:
        void nameFromTopology(const gmx::TopologyInformation &top);
        void chainFromTopology(const gmx::TopologyInformation &top);
        void hydrophobicityFromJson(const rapidjson::Document &doc);
        
        // getter methods:
        std::string name(const int id) const;
        std::string chain(const int id) const;
        real hydrophobicity(const int id) const;


    private:
    
        // container for residue properties:
        std::map<int, std::string> name_;
        std::map<int, std::string> chain_;
        std::map<int, real> hydrophobicity_;
};

#endif

