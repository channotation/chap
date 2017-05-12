#ifndef MOLECULAR_PATH_OBJ_EXPORTER_HPP
#define MOLECULAR_PATH_OBJ_EXPORTER_HPP

#include <gtest/gtest.h>

#include "io/obj_exporter.hpp"
#include "path-finding/molecular_path.hpp"


/*
 *
 */
class MolecularPathObjExporter : public ObjExporter
{
    friend class MolecularPathObjExporterTest;
    FRIEND_TEST(MolecularPathObjExporterTest, MolecularPathObjExporterOrthogonalVectorTest);
    FRIEND_TEST(MolecularPathObjExporterTest, MolecularPathObjExporterAxisRotationTest);


    public:
        
        // constructor:
        MolecularPathObjExporter();

        // interface for exporting:
        void operator()(char *filename, 
                        MolecularPath &molPath);

    private:

        //
        const real PI_ = std::acos(-1.0);

        // 
        inline int numPlanarVertices(real &d, real &r);

        // 
        gmx::RVec orthogonalVector(gmx::RVec vec);
        gmx::RVec rotateAboutAxis(gmx::RVec vec, gmx::RVec axis, real angle);
};


#endif

