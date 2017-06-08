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


/*!
 * This test checks that the radius provider has the correct behaviout in case
 * of broken input. In particular it checks that an exception is thrown if the
 * provided JSON document is broken, if it is missing an array containing the
 * van der Waals radius records, or if these records do not have the correct
 * entries with the correct types.
 */
TEST_F(VdwRadiusProviderTest, VdwRadiusProviderJsonTest)
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
    
    // create a JSON doc:
    rapidjson::Document jsonDocBrokenRecord;
    rapidjson::Document::AllocatorType& allocBrokenRecord = jsonDocBrokenRecord.GetAllocator();
    jsonDocBrokenRecord.SetObject();

    // create record with wrong type vdwr entry:
    rapidjson::Value brokenRecord1;
    brokenRecord1.SetObject();
    brokenRecord1.AddMember("atomname", "CA", allocBrokenRecord);
    brokenRecord1.AddMember("resname", "ARG", allocBrokenRecord);
    brokenRecord1.AddMember("vdwr", "1.0", allocBrokenRecord);

    // create record with wrong type atomname entry:
    rapidjson::Value brokenRecord2;
    brokenRecord2.SetObject();
    brokenRecord2.AddMember("atomname", 2.0, allocBrokenRecord);
    brokenRecord2.AddMember("resname", "ARG", allocBrokenRecord);
    brokenRecord2.AddMember("vdwr", 1.0, allocBrokenRecord);

    // create record with wrong type resname entry:
    rapidjson::Value brokenRecord3;
    brokenRecord3.SetObject();
    brokenRecord3.AddMember("atomname", "CA", allocBrokenRecord);
    brokenRecord3.AddMember("resname", -1.0, allocBrokenRecord);  
    brokenRecord3.AddMember("vdwr", 1.0, allocBrokenRecord);

    // create record with typo in vdwr entryi:
    rapidjson::Value brokenRecord4;
    brokenRecord4.SetObject();
    brokenRecord4.AddMember("atomname", "CA", allocBrokenRecord);
    brokenRecord4.AddMember("resname", "ARG", allocBrokenRecord);
    brokenRecord4.AddMember("vdwR", 1.0, allocBrokenRecord);

    // create record with typo in atomname entry:
    rapidjson::Value brokenRecord5;
    brokenRecord5.SetObject();
    brokenRecord5.AddMember("atmname", "CA", allocBrokenRecord);
    brokenRecord5.AddMember("resname", "ARG", allocBrokenRecord);
    brokenRecord5.AddMember("vdwr", 1.0, allocBrokenRecord);

    // create record with wrong type resname entry:
    rapidjson::Value brokenRecord6;
    brokenRecord6.SetObject();
    brokenRecord6.AddMember("atomname", "CA", allocBrokenRecord);
    brokenRecord6.AddMember("Resname", "ARG", allocBrokenRecord);  
    brokenRecord6.AddMember("vdwr", 1.0, allocBrokenRecord);

    // create array and add to document:
    rapidjson::Value vdwrarray(rapidjson::kArrayType);
    jsonDocBrokenRecord.AddMember("vdwradii", vdwrarray, allocBrokenRecord);

    // try creatings lookup tables from all these records:
    jsonDocBrokenRecord["vdwradii"].PushBack(brokenRecord1, allocBrokenRecord);
    ASSERT_THROW(rp.lookupTableFromJson(jsonDocBrokenRecord), std::runtime_error);
    jsonDocBrokenRecord["vdwradii"].Clear();

    jsonDocBrokenRecord["vdwradii"].PushBack(brokenRecord2, allocBrokenRecord);
    ASSERT_THROW(rp.lookupTableFromJson(jsonDocBrokenRecord), std::runtime_error);
    jsonDocBrokenRecord["vdwradii"].Clear();

    jsonDocBrokenRecord["vdwradii"].PushBack(brokenRecord3, allocBrokenRecord);
    ASSERT_THROW(rp.lookupTableFromJson(jsonDocBrokenRecord), std::runtime_error);
    jsonDocBrokenRecord["vdwradii"].Clear();

    jsonDocBrokenRecord["vdwradii"].PushBack(brokenRecord4, allocBrokenRecord);
    ASSERT_THROW(rp.lookupTableFromJson(jsonDocBrokenRecord), std::runtime_error);
    jsonDocBrokenRecord["vdwradii"].Clear();

    jsonDocBrokenRecord["vdwradii"].PushBack(brokenRecord5, allocBrokenRecord);
    ASSERT_THROW(rp.lookupTableFromJson(jsonDocBrokenRecord), std::runtime_error);
    jsonDocBrokenRecord["vdwradii"].Clear();

    jsonDocBrokenRecord["vdwradii"].PushBack(brokenRecord6, allocBrokenRecord);
    ASSERT_THROW(rp.lookupTableFromJson(jsonDocBrokenRecord), std::runtime_error);
    jsonDocBrokenRecord["vdwradii"].Clear();


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


/*
 *
 */
TEST_F(VdwRadiusProviderTest, VdwRadiusProviderLookupTest)
{

}

