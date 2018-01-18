#ifndef RESIDUE_INFORMATION_PROVIDER_HPP
#define RESIDUE_INFORMATION_PROVIDER_HPP

#include <map>
#include <string>
#include <vector>

#include <gromacs/topology/atoms.h>
#include <gromacs/topology/topology.h>
#include <gromacs/trajectoryanalysis/analysissettings.h>
#include <gromacs/utility/real.h>

#include "external/rapidjson/document.h"


/*!
 * Enum for hydrophobicity databases bundled with source code. 
 */
enum eHydrophobicityDatabase {eHydrophobicityDatabaseHessa2005,
                              eHydrophobicityDatabaseKyteDoolittle1982,
                              eHydrophobicityDatabaseMonera1995,
                              eHydrophobicityDatabaseMoon2011,
                              eHydrophobicityDatabaseWimleyWhite1996,
                              eHydrophobicityDatabaseZhu2016,
                              eHydrophobicityDatabaseMemprotMd,
                              eHydrophobicityDatabaseUser};


/*!
 * \brief Container class for (bio)chemical information about residues.
 *
 * This class provides information on residues by ID, such as residue name,
 * chain, or hydrophobicity.
 */
class ResidueInformationProvider
{
    public:

        // constructor and destructor:
        ResidueInformationProvider();

        // setter methods:
        void nameFromTopology(const gmx::TopologyInformation &top);
        void chainFromTopology(const gmx::TopologyInformation &top);
        void hydrophobicityFromJson(const rapidjson::Document &doc);
        void setDefaultHydrophobicity(const real hydrophobicity);
        
        // getter methods:
        std::vector<int> ids() const;
        std::string name(const int id) const;
        std::string chain(const int id) const;
        real hydrophobicity(const int id) const;


    private:
    
        // container for residue properties:
        std::map<int, std::string> name_;
        std::map<int, std::string> chain_;
        std::map<std::string, real> hydrophobicity_;

        // default properties:
        real defaultHydrophobicity_;
};

#endif

