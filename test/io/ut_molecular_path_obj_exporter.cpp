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


#include <vector>
#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "io/molecular_path_obj_exporter.hpp"


/*!
 * \brief Test fixture for the MolecularPathObjExporter.
 */
class MolecularPathObjExporterTest : public ::testing::Test
{
    public:
  
};



/*!
 * Test orthogonal vector construction.
 */
TEST_F(MolecularPathObjExporterTest, MolecularPathObjExporterOrthogonalVectorTest)
{
    // object to test on:
    MolecularPathObjExporter molPathExp;

    // create a set of test vectors:
    std::vector<gmx::RVec> vectors;
    vectors.push_back(gmx::RVec( 1.0,  0.0,  0.0));
    vectors.push_back(gmx::RVec( 0.0,  1.0,  0.0));
    vectors.push_back(gmx::RVec( 0.0,  0.0,  1.0));
    vectors.push_back(gmx::RVec(-1.0,  2.0, -3.5));

    // construct orthogonal vectors and check orthogonality:
    for(unsigned int i = 0; i < vectors.size(); i++)
    {
        gmx::RVec orthVec = molPathExp.orthogonalVector(vectors[i]);
        ASSERT_NEAR(0.0, iprod(vectors[i], orthVec), std::numeric_limits<real>::epsilon());
    }
}


/*!
 * Tests vector rotation around axis.
 */
TEST_F(MolecularPathObjExporterTest, MolecularPathObjExporterAxisRotationTest)
{
    const real PI = std::acos(-1.0);
    const real eps = std::numeric_limits<real>::epsilon();

    // object to test on:
    MolecularPathObjExporter molPathExp;

    // basis vectors:
    gmx::RVec vecX(1.0, 0.0, 0.0);
    gmx::RVec vecY(0.0, 1.0, 0.0);
    gmx::RVec vecZ(0.0, 0.0, 1.0);

    // rotate basis vectors around themselves:
    gmx::RVec rotX = molPathExp.rotateAboutAxis(vecX, vecX, PI);
    gmx::RVec rotY = molPathExp.rotateAboutAxis(vecY, vecY, PI);
    gmx::RVec rotZ = molPathExp.rotateAboutAxis(vecZ, vecZ, PI);
   
    // assert that they are unchanged:
    for(int i = 0; i < 3; i++)
    {
        ASSERT_NEAR(vecX[i], rotX[i], eps);
        ASSERT_NEAR(vecY[i], rotY[i], eps);
        ASSERT_NEAR(vecZ[i], rotZ[i], eps);
    }

    // rotate basis vectors by 90 degrees around one another:
    rotX = molPathExp.rotateAboutAxis(vecX, vecZ, PI/2.0);
    rotY = molPathExp.rotateAboutAxis(vecY, vecX, PI/2.0);
    rotZ = molPathExp.rotateAboutAxis(vecZ, vecY, PI/2.0);
 
    // assert that they now match up with the respective third axis:
    for(int i = 0; i < 3; i++)
    {
        ASSERT_NEAR(vecX[i], rotZ[i], eps);
        ASSERT_NEAR(vecY[i], rotX[i], eps);
        ASSERT_NEAR(vecZ[i], rotY[i], eps);
    }

    // rotate basis vectors by 180 degrees around one another:
    rotX = molPathExp.rotateAboutAxis(vecX, vecZ, PI);
    rotY = molPathExp.rotateAboutAxis(vecY, vecX, PI);
    rotZ = molPathExp.rotateAboutAxis(vecZ, vecY, PI);

    // assert that they are now mirrored:
    for(int i = 0; i < 3; i++)
    {
        ASSERT_NEAR(vecX[i], -rotX[i], eps);
        ASSERT_NEAR(vecY[i], -rotY[i], eps);
        ASSERT_NEAR(vecZ[i], -rotZ[i], eps);
    }

    // rotate a generic vector around each axis:
    gmx::RVec vec(1.0, -7.5, 3.1);
    rotX = molPathExp.rotateAboutAxis(vec, vecX, PI);
    rotY = molPathExp.rotateAboutAxis(vec, vecY, PI);
    rotZ = molPathExp.rotateAboutAxis(vec, vecZ, PI);
    
    // should be mirrored in other two components:
    ASSERT_NEAR( vec[XX], rotX[XX], 10*eps);
    ASSERT_NEAR(-vec[YY], rotX[YY], 10*eps);
    ASSERT_NEAR(-vec[ZZ], rotX[ZZ], 10*eps);
    ASSERT_NEAR(-vec[XX], rotY[XX], 10*eps);
    ASSERT_NEAR( vec[YY], rotY[YY], 10*eps);
    ASSERT_NEAR(-vec[ZZ], rotY[ZZ], 10*eps);
    ASSERT_NEAR(-vec[XX], rotZ[XX], 10*eps);
    ASSERT_NEAR(-vec[YY], rotZ[YY], 10*eps);
    ASSERT_NEAR( vec[ZZ], rotZ[ZZ], 10*eps);
}

