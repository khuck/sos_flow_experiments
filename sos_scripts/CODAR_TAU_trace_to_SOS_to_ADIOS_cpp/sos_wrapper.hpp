/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#pragma once

#include "utils.hpp"
#include "sos.h"
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

using json = nlohmann::json;

namespace extractor {

/* Class containing SOS connection info */
class sos {
    private:
        bool connected;
        json config;
        SOS_runtime* SOS;
    public:
        sos() : connected(false), SOS(nullptr) {
            PRINTSTACK
            read_config("./sos_config.json");
            //connect();
        };
        ~sos() {
            PRINTSTACK
            this->disconnect();
        };
        void connect() {
            PRINTSTACK
            SOS_init(&(this->SOS), SOS_ROLE_ANALYTICS, SOS_RECEIVES_NO_FEEDBACK, nullptr);
            if (SOS == nullptr) {
                fprintf(stderr, "ERROR: Could not initialize SOS.\n");
                exit (EXIT_FAILURE);
            }
            connected = true;
            srandom(SOS->my_guid);
        }
        void disconnect() {
            PRINTSTACK
            if (connected && this->SOS != nullptr) {
                SOS_finalize(this->SOS);
            }
        };
        void read_config(std::string filename) {
            PRINTSTACK
            // read a JSON file
            std::ifstream config_file(filename);
            config_file >> config;
            //std::cout << std::setw(2) << config << std::endl;
        }
        json& get_config(void) { return config; }
};

}; // end namespace extractor
