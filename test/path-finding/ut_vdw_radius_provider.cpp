#include <gtest/gtest.h>

#include "gromacs/topology/topology.h"

#include "io/json_doc_importer.hpp"
#include "path-finding/vdw_radius_provider.hpp"


/*
 *
 */
class VdwRadiusProviderTest : public ::testing::Test
{

};


/*
 *
 */
TEST_F(VdwRadiusProviderTest, VdwRadiusProviderSimpleTest)
{
    // read JSON document from file:    
    JsonDocImporter ImportJson;
    rapidjson::Document jsonDoc = ImportJson("../data/vdwradii/simple.json");

    // prepare radius provider: 
    VdwRadiusProvider rp;
    rp.lookupTableFromJson(jsonDoc);


   
    




    std::cout<<"rad = "<<rp.vdwRadiusForAtom("C", "??")<<std::endl;
}
