#ifndef OBJ_EXPORTER_HPP
#define OBJ_EXPORTER_HPP

#include <fstream>
#include <iostream>
#include <vector>

#include <gromacs/math/vec.h>
#include <gromacs/utility/real.h>   


/*!
 * Simple serialiser for writing data to a Wavefront OBJ geometry file. 
 * Currently, only writing of comments, vertices and faces is supported. Does
 * not perform error checking.
 */
class ObjExporter
{
    public:

        // interface for export:
        void write(char *filename,
                   std::vector<gmx::RVec> vertices,
                   std::vector<std::vector<int>> faces);


    private:

        // file handle:
        std::fstream obj_;

        // utilities for writing individual lines:
        inline void writeComment(char *comment);
        inline void writeVertex(gmx::RVec vertex);
        inline void writeFace(std::vector<int> face);
};

#endif

