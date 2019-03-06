/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#pragma once

#include "utils.hpp"
#include <nlohmann/json.hpp>
#include <string.h>

using json = nlohmann::json;

namespace extractor {

/* Class containing SOS connection info */
class sos {
    private:
        bool connected;
        json config;
        std::string hostname;
        int portnumber;
    public:
        sos(json& _config) : 
            connected(false),
            config(_config),
            hostname("localhost"),
            portnumber(22500)
        {
            PRINTSTACK
            connect();
        };
        ~sos() {
            PRINTSTACK
            disconnect();
        };
        void connect();
        void test_connection();
        void disconnect();
        bool check_for_frame(int frame);
        void sql_query(void);
        void cache_grab(void);
};

}; // end namespace extractor
