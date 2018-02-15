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


#ifndef VDW_RADIUS_PROVIDER_HPP
#define VDW_RADIUS_PROVIDER_HPP

#include <unordered_map>

#include <gtest/gtest.h>

#include <gromacs/trajectoryanalysis/analysissettings.h>
#include <gromacs/utility/arrayref.h>

#include "external/rapidjson/document.h"


/*!
 * \brief Abstract data type for van der Waals radius records.
 *
 * Bundles together atom name, residue name, and van der Waals radius and is 
 * using to build a lookup table in VdwRadiusProvider.
 */
struct VdwRadiusRecord
{
    std::string atmName_;
    std::string resName_;
    real vdwRad_;
};


/*!
 * \brief Enum for existing van-der-Waals radius databases distributed with the
 * source code.
 */
enum eVdwRadiusDatabase {eVdwRadiusDatabaseHoleAmberuni,
                         eVdwRadiusDatabaseHoleBondi,
                         eVdwRadiusDatabaseHoleHardcore,
                         eVdwRadiusDatabaseHoleSimple,
                         eVdwRadiusDatabaseHoleXplor,
                         eVdwRadiusDatabaseUser};


/*!
 * \brief Class for uniformly specifying van-der-Waals radii.
 *
 * This class is used to associate a van-der-Waals radius atoms in a given
 * topology. To this end, it internally maintains a lookup table of 
 * VdwRadiusRecord objects which are used to bundle together atom name, residue
 * name, and corresponding van-der-Waals radius. Such a lookup table can be
 * generated from a JSON document using lookupTableFromJson(). The 
 * dicumentation for this function explains the specifications of the JSON 
 * document.
 *
 * Once a lookup table has been created, vdwRadiiForTopology() can be used to
 * retrieve the van-der-Waals radii of atoms in that topology, where radii
 * are only returned for atoms whose ID is passed to vdwRadiiForTopology(). 
 * Internally, VdwRadiusProvider will look for matching records in its lookup
 * table according to the decision tree below:
 *
 *  \dot
 *  digraph example {
 *      node [shape=oval];
 *      a [ label="exact atom name match?"];
 *      b [ label="exact residue name match?"];
 *      c [ label="wildcard atom name match?"];
 *      d [ label="use value"];
 *      f [ label="wildcard residue name match?"];
 *      e [ label="use value"];
 *      h [ label="exact residue name match?"];
 *      i [ label="use value"];
 *      j [ label="wildcard residue name match?"];
 *      k [ label="use value"];
 *      l [ label="element name match?"];
 *      m [ label="exact residue name match?"];
 *      n [ label="use value"];
 *      o [ label="wildcard residue name match?"];
 *      p [ label="use value"];
 *      q [ label="default radius available?"];
 *      r [ label="use value"];
 *      s [ label="throw exception"];
 *      a -> b [ arrowhead="open", label = "yes"];
 *      a -> c [ arrowhead="open", label = "no"];
 *      b -> d [ arrowhead="open", label = "yes"];
 *      b -> f [ arrowhead="open", label = "no"];
 *      f -> e [ arrowhead="open", label = "yes"];
 *      f -> c [ arrowhead="open", label = "no"];
 *      c -> h [ arrowhead="open", label = "yes"];
 *      h -> i [ arrowhead="open", label = "yes"];
 *      h -> j [ arrowhead="open", label = "no"];
 *      j -> k [ arrowhead="open", label = "yes"];
 *      c -> l [ arrowhead="open", label = "no"];
 *      j -> l [ arrowhead="open", label = "no"];
 *      l -> m [ arrowhead="open", label = "yes"];
 *      m -> n [ arrowhead="open", label = "yes"];
 *      m -> o [ arrowhead="open", label = "no"];
 *      o -> p [ arrowhead="open", label = "yes"];
 *      o -> q [ arrowhead="open", label = "no"];
 *      l -> q [ arrowhead="open", label = "no"];
 *      q -> r [ arrowhead="open", label = "yes"];
 *      q -> s [ arrowhead="open", label = "no"];
 *  }
 * \enddot
 *
 * The VdwRadiusProvider will always try to find an exactly matching atom name
 * first, before trying to find a matching wildcard atom name or a matching
 * element name. A wildcard atom name is one where some characters in the atom
 * name have been replaced by questions marks, e.g. 'C???' would be a wildcard
 * for an atom whose name begins with 'C' and may contain up to three further
 * arbitrary characters. If, for example, a topology contains an atom named 
 * 'CA' it would not be a match for record named 'C', but would be a wildcard
 * match for a record named 'C???' (it would also be a match for a record named
 * 'CA'). If VdwRadiusProvider fails to find a record with an atom name that
 * matches that of the trial atom either exactly or with wildcards, it will try
 * to find a record that matches the atom's (capitalised) element name, i.e.
 * if there are no records named 'CE2', VdwRadiusProvider will try to find a 
 * record named 'C' even if no wildcard like 'C???' is in the lookup table.
 *
 * VdwRadiusProvider will internally collect all records with matching atom
 * name (or matching wildcard atom name or matching element name) and will then
 * try to find a record with matching residue name amongst these. Again, an
 * exact name match is attempted first before trying a wildcard match. In 
 * contrast to the atom name case, resdiue name wildcards always need to 
 * contain exactly three questionmarks (i.e. '???'), partial wildcards (like 
 * 'G??' will be overlooked. If neither an exact or a wildcard residue name 
 * match is found, VdwRadiusProvider will continue its lookup with trying to
 * match wildcard atom names (if it failed to find a matching residue name for
 * a set of matching atom names) and with trying to find matching element names
 * (if it failed to find a matching residue name for a set of matching wildcard
 * atom names).
 *
 * If no matching record could be found after trying to element names, the last
 * resort is to use default radius that can be set with setDefaultVdwRadius()
 * prior to calling vdwRadiiForTopology(). If not default radius has been set,
 * the lookup will throw an exception.
 */
class VdwRadiusProvider
{
    friend class VdwRadiusProviderTest;
    FRIEND_TEST(VdwRadiusProviderTest, VdwRadiusProviderJsonTest);
    FRIEND_TEST(VdwRadiusProviderTest, VdwRadiusProviderLookupTest);

    public:

        // constructor and destructor:
        VdwRadiusProvider();
        ~VdwRadiusProvider();

        // public interface for setting a default radius:
        void setDefaultVdwRadius(real defRad);

        // public interface for generating vdW radius lookup table:
        void lookupTableFromJson(rapidjson::Document &jsonDoc);

        // public interface for obtaining vdwRadii for given topology:
        std::unordered_map<int, real> vdwRadiiForTopology(
            const gmx::TopologyInformation &top,
            std::vector<int> mappedIds);

    private:

        // default vdW radius to be used if no match in lookup table:
        real defRad_;

        // lookup table for vdW radii:
        std::vector<VdwRadiusRecord> vdwRadiusLookupTable_;

        // function to perform sanity checks on lookup table:
        void validateLookupTable();

        // function for associating a vdW radius with an atom and residue name:
        real vdwRadiusForAtom(std::string atmName, 
                              std::string resName,
                              std::string elemSym);
        std::vector<VdwRadiusRecord> matchAtmName(std::string atmName);
        std::vector<VdwRadiusRecord> matchPartAtmName(std::string atmName);
        std::vector<VdwRadiusRecord>::const_iterator matchResName(
            std::string resName,
            const std::vector<VdwRadiusRecord> &records);

        // function to validate and return default radius:
        inline real returnDefaultRadius(std::string atmName, std::string resName);
};

#endif

