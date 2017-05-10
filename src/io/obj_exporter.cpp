#include "io/obj_exporter.hpp"


/*!
 * Public interface for OBJ export. Takes a file name and lists of vertices
 * and faces, which are then written to an OBJ format file. 
 */
void
ObjExporter::write(char *filename,
                   std::vector<gmx::RVec> vertices,
                   std::vector<std::vector<int>> faces)
{
    // open file steam:
    obj_.open(filename, std::fstream::out);

    // writer header comment:
    writeComment("produced by CHAP");
        
    // write vertices:
    obj_<<std::endl;
    for(unsigned int i = 0; i < vertices.size(); i++)
    {
        writeVertex(vertices[i]);
    }

    // write faces: 
    obj_<<std::endl;
    for(unsigned int i = 0; i < faces.size(); i++)
    {
        writeFace(faces[i]); 
    }

    // close file stream:
    obj_.close();
}


/*!
 * Writes a comment line to an OBJ file.
 */
void
ObjExporter::writeComment(char *comment)
{
    obj_ <<"# "<<comment<<std::endl;
}


/*!
 * Writes a vertex entry to an OBJ file.
 */
void
ObjExporter::writeVertex(gmx::RVec vertex)
{
    obj_<<"v "<<vertex[0]<<" "
              <<vertex[1]<<" "
              <<vertex[2]<<std::endl;
}


/*!
 * Writes a face entry to OBJ file.
 */
void
ObjExporter::writeFace(std::vector<int> face)
{
    obj_<<"f ";

    for(unsigned int i = 0; i < face.size(); i++)
    {
        obj_<<face[i]<<"  ";
    }

    obj_<<std::endl;
}

