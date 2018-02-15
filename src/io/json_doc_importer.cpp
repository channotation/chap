#include <fstream>
#include <exception>
#include <stdexcept>
#include <sstream>
#include <iostream>

#include "io/json_doc_importer.hpp"


/*!
 * Returns a JSON document corresponding to the given JSON file.
 */
rapidjson::Document
JsonDocImporter::operator()(std::string fileName)
{
    // read entire file into memory:
    std::stringstream ss;
    std::ifstream file;
    file.open(fileName.c_str());

    // make sure file could be opened:
    if( !file.is_open() )
    {
        throw std::runtime_error("ERROR: Could not open file " + fileName + ".");
    }

    // read file into string stream:
    ss << file.rdbuf();
    file.close();

    // create JSON document from file:
    rapidjson::Document json;
    json.Parse<0>(ss.str().c_str());

    // check validity of JSON object:
    if( json.IsObject() == false )
    {
        throw std::invalid_argument("Invalid JSON object.");
    }

    // return JSON document:
    return json;
}

