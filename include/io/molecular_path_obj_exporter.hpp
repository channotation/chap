#ifndef MOLECULAR_PATH_OBJ_EXPORTER_HPP
#define MOLECULAR_PATH_OBJ_EXPORTER_HPP

#include "io/obj_exporter.hpp"
#include "path-finding/molecular_path.hpp"


/*
 *
 */
class MolecularPathObjExporter : public ObjExporter
{
    public:
        
        // constructor:
        MolecularPathObjExporter();

        // interface for exporting:
        void operator()(char *filename, 
                        MolecularPath molPath);

    private:

        //
        const real PI_ = std::acos(-1.0);

        // 
        inline int numPlanarVertices(real &d, real &r);

};


#endif

