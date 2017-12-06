#ifndef WAVEFRONT_MTL_IO_HPP
#define WAVEFRONT_MTL_IO_HPP

#include <fstream>
#include <string>
#include <vector>

#include <gromacs/math/vec.h>
#include <gromacs/utility/real.h>


/*
 *
 */
class WavefrontMtlMaterial
{
    friend class WavefrontMtlExporter;

    public:

        // setter functions for material properties:
        void setAmbientColour(gmx::RVec col);
        void setDiffuseColour(gmx::RVec col);
        void setSpecularColour(gmx::RVec col);

    private:

        // material properties:
        std::string name_;
        gmx::RVec ambientColour_;
        gmx::RVec diffuseColour_;
        gmx::RVec specularColour_;

};


/*
 *
 */
class WavefrontMtlObject
{
    friend class WavefrontMtlExporter;

    public:

        // interface for adding materials:
        void addMaterial(WavefrontMtlMaterial material);

    private:

        // collection of materials:
        std::vector<WavefrontMtlMaterial> materials_;
};


/*
 *
 */
class WavefrontMtlExporter
{
    public:

        // public interfacr for writing material definitions to file:
        void write(std::string fileName, WavefrontMtlObject object);

    private:

        // file handle:
        std::fstream file_;
    
        // utilities for writing individual lines:
        void writeMaterialName(std::string name);
        void writeAmbientColour(const gmx::RVec &col);
        void writeDiffuseColour(const gmx::RVec &col);
        void writeSpecularColour(const gmx::RVec &col);
};

#endif

