/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#pragma once

#include "utils.hpp"
#include <nlohmann/json.hpp>
#include "adios2.h"

using json = nlohmann::json;

namespace extractor {

/* Class containing ADIOS archive info */
class adios {
    private:
        bool opened;
        json config;
        adios2::ADIOS ad;
        adios2::IO bpIO;
        adios2::Engine bpWriter;
    public:
        adios(json& _config) : 
            opened(false),
            config(_config)
        {
            PRINTSTACK
            initialize();
            open();
        };
        ~adios() {
            PRINTSTACK
            close();
        };
        void initialize();
        void open();
        void close();
};

}; // end namespace extractor
