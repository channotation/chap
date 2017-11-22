#ifndef MOLECULAR_PATH_OBJ_EXPORTER_HPP
#define MOLECULAR_PATH_OBJ_EXPORTER_HPP

#include <map>
#include <string>

#include <gtest/gtest.h>

#include "io/wavefront_obj_io.hpp"
#include "path-finding/molecular_path.hpp"


/*
 *
 */
class RegularVertexGrid
{
    public:

        RegularVertexGrid(
                std::vector<real> s,
                std::vector<real> phi);

        void addVertex(size_t i, size_t j, gmx::RVec vertex);

        void interpolateMissing();
    
        std::vector<gmx::RVec> vertices();
        std::vector<gmx::RVec> normals();
        std::vector<WavefrontObjFace> faces();

    private:

        int numX;
        int numY;

        const std::vector<real> phi_;
        const std::vector<real> s_;

        std::map<std::pair<size_t, size_t>, gmx::RVec> vertices_;
//        std::vector<gmx::RVec> vertices_;
        
};


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

