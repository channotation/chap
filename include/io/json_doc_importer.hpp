#ifndef JSON_DOC_IMPORTER_HPP
#define JSON_DOC_IMPORTER_HPP

#include <string>

#include "external/rapidjson/document.h"


/*
 *
 */
class JsonDocImporter
{
    public:

        // constructor and destructor:
        JsonDocImporter(){};
        ~JsonDocImporter(){};

        // define operator for file reading:
        rapidjson::Document operator()(std::string fileName);
};

#endif

