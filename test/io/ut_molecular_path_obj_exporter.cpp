#include <vector>
#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "io/molecular_path_obj_exporter.hpp"


/*
 *
 */
class MolecularPathObjExporterTest : public ::testing::Test
{
    public:
  
};



/*
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

    // construct orthgonal vectors and check orthogonality:
    for(unsigned int i = 0; i < vectors.size(); i++)
    {
        gmx::RVec orthVec = molPathExp.orthogonalVector(vectors[i]);
        ASSERT_NEAR(0.0, iprod(vectors[i], orthVec), std::numeric_limits<real>::epsilon());
    }
}


/*
 *
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
    ASSERT_NEAR( vec[0], rotX[0], 10*eps);
    ASSERT_NEAR(-vec[1], rotX[1], 10*eps);
    ASSERT_NEAR(-vec[2], rotX[2], 10*eps);
    ASSERT_NEAR(-vec[0], rotY[0], 10*eps);
    ASSERT_NEAR( vec[1], rotY[1], 10*eps);
    ASSERT_NEAR(-vec[2], rotY[2], 10*eps);
    ASSERT_NEAR(-vec[0], rotZ[0], 10*eps);
    ASSERT_NEAR(-vec[1], rotZ[1], 10*eps);
    ASSERT_NEAR( vec[2], rotZ[2], 10*eps);
}

