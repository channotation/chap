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


#include <limits>

#include <gtest/gtest.h>

#include "gromacs/topology/topology.h"

#include "io/json_doc_importer.hpp"
#include "path-finding/vdw_radius_provider.hpp"


/*!
 * \brief Test fixture for VdwRadiusProvider. 
 *
 * Only used to group test cases.
 */
class VdwRadiusProviderTest : public ::testing::Test
{

};


/*!
 * This test checks that the radius provider has the correct behaviour in case
 * of broken input. In particular it checks that an exception is thrown if the
 * provided JSON document is broken, if it is missing an array containing the
 * van der Waals radius records, or if these records do not have the correct
 * entries with the correct types. It also checks that an exception is thworn 
 * if a lookup table with duplicate entries is created.
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


    // DUPLICATE RECORD CASE
    //-------------------------------------------------------------------------

    rapidjson::Document jsonDocDuplicate;
    rapidjson::Document::AllocatorType& allocDuplicate = jsonDocDuplicate.GetAllocator();

    // create correct, but duplicate entries:
    rapidjson::Value duplicateRecord1;
    duplicateRecord1.SetObject();
    duplicateRecord1.AddMember("atomname", "C", allocDuplicate);
    duplicateRecord1.AddMember("resname", "ARG", allocDuplicate);
    duplicateRecord1.AddMember("vdwr", 1.0, allocDuplicate);

    rapidjson::Value duplicateRecord2;
    duplicateRecord2.SetObject();
    duplicateRecord2.AddMember("atomname", "CA", allocDuplicate);
    duplicateRecord2.AddMember("resname", "ARG", allocDuplicate);
    duplicateRecord2.AddMember("vdwr", 1.0, allocDuplicate);

    rapidjson::Value duplicateRecord3;
    duplicateRecord3.SetObject();
    duplicateRecord3.AddMember("atomname", "N", allocDuplicate);
    duplicateRecord3.AddMember("resname", "GLY", allocDuplicate);
    duplicateRecord3.AddMember("vdwr", 1.0, allocDuplicate);

    rapidjson::Value duplicateRecord4;
    duplicateRecord4.SetObject();
    duplicateRecord4.AddMember("atomname", "C", allocDuplicate);
    duplicateRecord4.AddMember("resname", "ARG", allocDuplicate);
    duplicateRecord4.AddMember("vdwr", 1.0, allocDuplicate);

    rapidjson::Value duplicateRecord5;
    duplicateRecord5.SetObject();
    duplicateRecord5.AddMember("atomname", "C", allocDuplicate);
    duplicateRecord5.AddMember("resname", "GLY", allocDuplicate);
    duplicateRecord5.AddMember("vdwr", 1.0, allocDuplicate);
    
    // add object to an array:
    rapidjson::Value vdwradiiDuplicate(rapidjson::kArrayType);
    vdwradiiDuplicate.PushBack(duplicateRecord1, allocDuplicate);
    vdwradiiDuplicate.PushBack(duplicateRecord2, allocDuplicate);
    vdwradiiDuplicate.PushBack(duplicateRecord3, allocDuplicate);
    vdwradiiDuplicate.PushBack(duplicateRecord4, allocDuplicate);
    vdwradiiDuplicate.PushBack(duplicateRecord5, allocDuplicate);

    // add array to document:
    jsonDocDuplicate.SetObject();
    jsonDocDuplicate.AddMember("vdwradii", vdwradiiDuplicate, allocDuplicate);

    // proper document should not throw:
    ASSERT_THROW(rp.lookupTableFromJson(jsonDocDuplicate), std::runtime_error);


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


/*!
 * This test case checks that the van der Waals radius provider finds the 
 * correct radius for any possible type of combination of input atom name, 
 * residue name, and element name. It starts by building a lookup table for
 * a few atom and residue name combinations and then checks that the correct
 * radius is returned in the following cases:
 *
 *  - match for atom name and exact residue name
 *  - match for atom name and generic residue name
 *  - match for element name and exact residue name
 *  - match for element name and generic residue name
 *  - no match for atom or residue name, but default radius set
 *
 * The test also makes sure that an exception is thrown in the following cases:
 *
 *  - match for atom name, but no match for residue name and no default set
 *  - match for element name, but noch match for residue name and no default set
 *  - no match for atom or alement name, and no default set
 */
TEST_F(VdwRadiusProviderTest, VdwRadiusProviderLookupTest)
{
    // floating point comparison threshold:
    real eps = std::numeric_limits<real>::epsilon();

    // test value for radii:
    std::vector<std::string> atomNames = {"P", "C", "C", "O", "CL", "N", "FE", "O???", "O???"};
    std::vector<std::string> resNames = {"GLY", "ARG", "???", "ARG", "CL", "???", "???", "LYS", "???"};
    std::vector<real> vdwRadii = {3.1, 1.1, 2.2, 3.3, 4.4, 6.6, 1.7, 2.8, 2.9};

    // create JSON document and allocator:
    rapidjson::Document radii;
    radii.SetObject();
    rapidjson::Document::AllocatorType& alloc = radii.GetAllocator();
   
    // add array of radius value to document:
    rapidjson::Value vdwrarray(rapidjson::kArrayType);
    for(size_t i = 0; i < vdwRadii.size(); i++)
    {
        rapidjson::Value record;
        record.SetObject();
        rapidjson::Value atomName;
        rapidjson::Value resName; 
        atomName.SetString(atomNames[i].c_str(), atomNames[i].size(), alloc);
        resName.SetString(resNames[i].c_str(), resNames[i].size(), alloc);
        record.AddMember("atomname", atomName, alloc);
        record.AddMember("resname", resName, alloc);
        record.AddMember("vdwr", vdwRadii[i], alloc);
        vdwrarray.PushBack(record, alloc);
    }
    radii.AddMember("vdwradii", vdwrarray, alloc);

    // create radius provider and build lookup table:
    VdwRadiusProvider rp;
    rp.lookupTableFromJson(radii);

    // check exact match for atom and residu name:
    ASSERT_NEAR(1.1, rp.vdwRadiusForAtom("C", "ARG", "C"), eps);
    ASSERT_NEAR(4.4, rp.vdwRadiusForAtom("CL", "CL", "Cl"), eps);

    // check exact match for atom name only:
    ASSERT_NEAR(2.2, rp.vdwRadiusForAtom("C", "TYR", "C"), eps);
    ASSERT_NEAR(6.6, rp.vdwRadiusForAtom("N", "ALA", "N"), eps);

    // check exact match for atom name only with no default residue:
    ASSERT_THROW(rp.vdwRadiusForAtom("P", "TYR", "P"), std::runtime_error);
    ASSERT_THROW(rp.vdwRadiusForAtom("NA", "NA", "Na"), std::runtime_error);

    // check case with partial atom name match: 
    ASSERT_NEAR(2.8, rp.vdwRadiusForAtom("O2", "LYS", "O"), eps);
    ASSERT_NEAR(2.9, rp.vdwRadiusForAtom("O2", "ARG", "O"), eps);

    // check case with element name and exact residue name match: 
    ASSERT_NEAR(1.1, rp.vdwRadiusForAtom("CA", "ARG", "C"), eps);
    ASSERT_NEAR(1.1, rp.vdwRadiusForAtom("CB", "ARG", "C"), eps);

    // check case with element name and generic residue match:
    ASSERT_NEAR(2.2, rp.vdwRadiusForAtom("CA", "TYR", "C"), eps);
    ASSERT_NEAR(1.7, rp.vdwRadiusForAtom("FE2", "TYR", "Fe"), eps);

    // check case with element name match but no residue name match:
    ASSERT_THROW(rp.vdwRadiusForAtom("P2", "TYR", "P"), std::runtime_error);

    // check no match with no default radius set:
    ASSERT_THROW(rp.vdwRadiusForAtom("E2", "ARG", "H"), std::runtime_error);
    ASSERT_THROW(rp.vdwRadiusForAtom("MW", "SOL", ""), std::runtime_error);
    ASSERT_THROW(rp.vdwRadiusForAtom("LP", "SOL", ""), std::runtime_error);

    // check no match with default radius set:
    rp.setDefaultVdwRadius(5.5);
    ASSERT_NEAR(5.5, rp.vdwRadiusForAtom("E2", "ARG", "H"), eps);
}

