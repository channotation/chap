#ifndef WAVEFRONT_OBJ_IO_HPP
#define WAVEFRONT_OBJ_IO_HPP

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gromacs/math/vec.h>
#include <gromacs/utility/real.h>   


/*!
 * Effectively an abstract data type representing groups OBJ files. Consists
 * of a group name and a vector of integer vectors representing the list of 
 * faces.
 */
class WavefrontObjGroup
{
    public:

        // constructor:
        WavefrontObjGroup(std::string name,
                          std::vector<std::vector<int>> faces);
        WavefrontObjGroup(const WavefrontObjGroup &other);

        // data:
        std::string groupname_;
        std::vector<std::vector<int>> faces_;
};


/*!
 * An abstract data type representing a complete Wavefront OBJ object (i.e. an
 * entire file). Consists of a vector of vertex positions and a vector of
 * WavefrongObjGroup objects.
 */
class WavefrontObjObject
{
    public:

        // constructor:
        WavefrontObjObject(std::string name);

        // functions to add data:
        void addVertices(std::vector<gmx::RVec> vertices);
        void addGroup(std::string name, std::vector<std::vector<int>> faces);

        // functions to manipulate data:
        void scale(real fac);
        void shift(gmx::RVec shift);

        // functions to query data:
        gmx::RVec calculateCog();

        // data:
        std::string name_;
        std::vector<gmx::RVec> vertices_;
        std::vector<WavefrontObjGroup> groups_;
};


/*!
 * Simple serialiser for writing data to a Wavefront OBJ geometry file. 
 * Currently, only writing of comments, vertices and faces is supported. Faces
 * may be grouped together. Does not perform error checking.
 */
class WavefrontObjExporter
{
    public:

        // interface for export:
        void write(std::string fileName,
                   std::vector<gmx::RVec> vertices,
                   std::vector<std::vector<int>> faces);

        void write(std::string fileName,
                   WavefrontObjObject object);


    private:

        // file handle:
        std::fstream obj_;

        // utilities for writing individual lines:
        inline void writeComment(std::string comment);
        inline void writeGroup(std::string group);
        inline void writeVertex(gmx::RVec vertex);
        inline void writeFace(std::vector<int> face);
};

#endif

