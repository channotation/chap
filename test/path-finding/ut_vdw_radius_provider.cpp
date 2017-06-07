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
    // create radius provider:
    VdwRadiusProvider rp;

    // BROKEN DOCUMENT CASE
    //-------------------------------------------------------------------------

    // using invalid JSON object should fail:
    rapidjson::Document jsonDocMissingData;
    rapidjson::Document::AllocatorType& allocMissingData = jsonDocMissingData.GetAllocator();
    ASSERT_THROW(rp.lookupTableFromJson(jsonDocMissingData), std::runtime_error);

    // using empty JSON object should fail:
    jsonDocMissingData.SetObject();
    ASSERT_THROW(rp.lookupTableFromJson(jsonDocMissingData), std::runtime_error);

    // using JSON with wrong type vdwradii should fail:
    jsonDocMissingData.AddMember("vdwradii", "wrong type", allocMissingData);
    ASSERT_THROW(rp.lookupTableFromJson(jsonDocMissingData), std::runtime_error);


    // BROKEN RADIUS RECORDS CASE
    //-------------------------------------------------------------------------
    
    // TODO: this!

    



    // WORKING CASE
    //-------------------------------------------------------------------------

    rapidjson::Document jsonDocNoThrow;
    rapidjson::Document::AllocatorType& allocNoThrow = jsonDocNoThrow.GetAllocator();

    // create an object with correct entries:
    rapidjson::Value record;
    record.SetObject();
    record.AddMember("atomname", "CA", allocNoThrow);
    record.AddMember("resname", "ARG", allocNoThrow);
    record.AddMember("vdwr", 1.0, allocNoThrow);
    
    // add object to an array:
    rapidjson::Value vdwradii(rapidjson::kArrayType);
    vdwradii.PushBack(record, allocNoThrow);

    // add array to document:
    jsonDocNoThrow.SetObject();
    jsonDocNoThrow.AddMember("vdwradii", vdwradii, allocNoThrow);

    // proper document should not throw:
    ASSERT_NO_THROW(rp.lookupTableFromJson(jsonDocNoThrow));
}
