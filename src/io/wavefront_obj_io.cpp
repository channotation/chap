#include "io/wavefront_obj_io.hpp"


/*!
 * Construct a new group with given name and list of faces.
 */
WavefrontObjGroup::WavefrontObjGroup(std::string name, 
                                     std::vector<std::vector<int>> faces)
    : groupname_(name)
    , faces_(faces)
{
    
}


/*!
 * Copy constructor for allowing use of STL containers.
 */
WavefrontObjGroup::WavefrontObjGroup(const WavefrontObjGroup &other)
    : groupname_(other.groupname_)
    , faces_(other.faces_)
{
    
}


/*!
 * Constructs and empty OBJ object with a given name.
 */
WavefrontObjObject::WavefrontObjObject(std::string name)
    : name_(name)
{
    
}


/*!
 * Add new vertices to the OBJ object. These are appended to the existing list
 * of vertices and there is no check for redundancy.
 */
void
WavefrontObjObject::addVertices(std::vector<gmx::RVec> vertices)
{
    vertices_.insert(vertices_.end(), vertices.begin(), vertices.end());
}


/*!
 * Adds a new named group to the OBJ object.
 */
void
WavefrontObjObject::addGroup(std::string name, std::vector<std::vector<int>> faces)
{
    groups_.push_back( WavefrontObjGroup(name, faces) );
}


/*
 *
 */
void
WavefrontObjObject::scale(real fac)
{
    // find centre of geometry:
    gmx::RVec cog = calculateCog();
    gmx::RVec negCog;
    negCog[0] = -cog[0];
    negCog[1] = -cog[1];
    negCog[2] = -cog[2];

    // shift vertices to be centred around origin:
    this -> shift(negCog);

    // scale all position vectors:
    std::vector<gmx::RVec>::iterator it;
    for(it = vertices_.begin(); it != vertices_.end(); it++)
    {
        (*it)[0] *= fac;
        (*it)[1] *= fac;
        (*it)[2] *= fac;
    }

    // shift vertices back to original centre of geometry:
    this -> shift(cog);
}


/*
 *
 */
void WavefrontObjObject::shift(gmx::RVec shift)
{
    std::vector<gmx::RVec>::iterator it;
    for(it = vertices_.begin(); it != vertices_.end(); it++)
    {
        (*it)[0] += shift[0];
        (*it)[1] += shift[1];
        (*it)[2] += shift[2];
    }    
}


/*
 *
 */
gmx::RVec
WavefrontObjObject::calculateCog()
{
    gmx::RVec cog(0.0, 0.0, 0.0);

    std::vector<gmx::RVec>::iterator it;
    for(it = vertices_.begin(); it != vertices_.end(); it++)
    {
        cog[0] += (*it)[0];
        cog[1] += (*it)[1];
        cog[2] += (*it)[2];
    }

    cog[0] /= vertices_.size();
    cog[1] /= vertices_.size();
    cog[2] /= vertices_.size();


    std::cout<<"cog = "<<cog[0]<<"  "<<cog[1]<<"  "<<cog[2]<<std::endl;

    return cog;
}


/*!
 * Writes an OBJ object to a file of the given name. 
 */
void
WavefrontObjExporter::write(char *filename,
                            WavefrontObjObject object)
{
    // open file stream:
    obj_.open(filename, std::fstream::out);

    // writer header comment:
    writeComment("produced by CHAP");
        
    // write vertices:
    obj_<<std::endl;
    for(unsigned int i = 0; i < object.vertices_.size(); i++)
    {
        writeVertex(object.vertices_[i]);
    }

    // write groups:
    std::vector<WavefrontObjGroup>::iterator it;
    for(it = object.groups_.begin(); it != object.groups_.end(); it++)
    {
        // write group name:
        writeGroup(it -> groupname_);

        // write faces in this group:
        for(int j = 0; j < it -> faces_.size(); j++)
        {
            writeFace(it -> faces_[j]);
        }
    }

    // close file stream:
    obj_.close();
}


/*!
 * Writes a comment line to an OBJ file.
 */
void
WavefrontObjExporter::writeComment(char *comment)
{
    obj_ <<"# "<<comment<<std::endl;
}

/*!
 * Writes a group line to an OBJ file.
 */
void
WavefrontObjExporter::writeGroup(std::string group)
{
    obj_ <<std::endl;
    obj_ <<"g "<<group<<std::endl;
}

/*!
 * Writes a vertex entry to an OBJ file.
 */
void
WavefrontObjExporter::writeVertex(gmx::RVec vertex)
{
    obj_<<"v "<<vertex[0]<<" "
              <<vertex[1]<<" "
              <<vertex[2]<<std::endl;
}

/*!
 * Writes a face entry to OBJ file.
 */
void
WavefrontObjExporter::writeFace(std::vector<int> face)
{
    obj_<<"f ";

    for(unsigned int i = 0; i < face.size(); i++)
    {
        obj_<<face[i]<<"  ";
    }

    obj_<<std::endl;
}

