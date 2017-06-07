#ifndef VDW_RADIUS_PROVIDER_HPP
#define VDW_RADIUS_PROVIDER_HPP

#include <unordered_map>

#include <gtest/gtest.h>

#include <gromacs/trajectoryanalysis/analysissettings.h>

#include "rapidjson/document.h"



/*
 *
 */
typedef struct VdwRadiusRecord
{
    std::string atmName_;
    std::string resName_;
    real vdwRad_;
};


/*
 *
 */
class VdwRadiusProvider
{
    friend class VdwRadiusProviderTest;
    FRIEND_TEST(VdwRadiusProviderTest, VdwRadiusProviderSimpleTest);

    public:

        // constructor and destructor:
        VdwRadiusProvider();
        ~VdwRadiusProvider();

        // public interface for setting a default radius:
        void setDefaultVdwRadius(real defRad);

        // public interface for generating vdW radius lookup table:
        void lookupTableFromJson(rapidjson::Document &jsonDoc);

        // public interface for obtaining vdwRadii for given topology:
        std::unordered_map<int, real> vdwRadiiForTopology(const gmx::TopologyInformation &top);

    private:

        // default vdW radius to be used if no match in lookup table:
        real defRad_;

        // lookup table for vdW radii:
        std::vector<VdwRadiusRecord> vdwRadiusLookupTable_;

        // function for associating a vdW radius with an atom and residue name:
        real vdwRadiusForAtom(std::string atmName, std::string resName);

        // function to validate and return default radius:
        inline real returnDefaultRadius(std::string atmName, std::string resName);
};

#endif

