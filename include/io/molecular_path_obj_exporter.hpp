#ifndef MOLECULAR_PATH_OBJ_EXPORTER_HPP
#define MOLECULAR_PATH_OBJ_EXPORTER_HPP

#include <string>

#include <gtest/gtest.h>

#include "path-finding/molecular_path.hpp"


/*
 *
 */
class MolecularPathObjExporter
{
    friend class MolecularPathObjExporterTest;
    FRIEND_TEST(MolecularPathObjExporterTest, MolecularPathObjExporterOrthogonalVectorTest);
    FRIEND_TEST(MolecularPathObjExporterTest, MolecularPathObjExporterAxisRotationTest);


    public:
        
        // constructor:
        MolecularPathObjExporter();

        // interface for exporting:
        void operator()(std::string fileName, 
                        MolecularPath &molPath);

    private:

        //
        const real PI_ = std::acos(-1.0);

        // 
        inline int numPlanarVertices(real &d, real &r);

        // 
        inline std::pair<std::vector<gmx::RVec>, std::vector<gmx::RVec>> vertexRing(
                gmx::RVec base,
                gmx::RVec tangent,
                gmx::RVec normal,
                real radius,
                real angleIncrement,
                size_t nIncrements);

        // 
        gmx::RVec orthogonalVector(gmx::RVec vec);
        gmx::RVec rotateAboutAxis(gmx::RVec vec, gmx::RVec axis, real angle);
};


#endif

