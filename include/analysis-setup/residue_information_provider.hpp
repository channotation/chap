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

