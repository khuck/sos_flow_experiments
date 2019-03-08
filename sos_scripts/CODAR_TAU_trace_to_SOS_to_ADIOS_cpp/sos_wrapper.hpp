/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#pragma once

#include "utils.hpp"
#include "adios_wrapper.hpp"
#include <nlohmann/json.hpp>
#include <string.h>
#include <unordered_map>
#include <vector>

using json = nlohmann::json;

namespace extractor {

/* Class containing SOS connection info */
class sos {
    private:
        bool connected;
        json config;
        std::string hostname;
        int portnumber;
        std::unordered_map<std::string, int> column_map;
        int prog_name_index;
        int comm_rank_index;
        int value_name_index;
        int value_index;
        int frame_index;
        int total_valid;
        int time_index;
        std::unordered_map<std::string, int> prog_names;
        std::unordered_map<std::string, int> comm_ranks;
        std::unordered_map<std::string, int> value_names;
        std::unordered_map<std::string, int> threads;
        std::unordered_map<std::string, int> metadata_keys;
        std::unordered_map<std::string, int> groups;
        std::unordered_map<std::string, int> timers;
        std::unordered_map<std::string, int> event_types;
        static std::vector<std::string> split(const std::string& s, char delimiter);
    public:
        sos(json& _config) : 
            connected(false),
            config(_config),
            hostname("localhost"),
            portnumber(22500),
            prog_name_index(-1),
            comm_rank_index(-1),
            value_name_index(-1),
            value_index(-1),
            frame_index(-1),
            time_index(-1)
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
        void build_column_map(void);
        void write_metadata(int frame, adios& my_adios);
        void write_timer_data(int frame, adios& my_adios);
        void sql_query(void);
        void cache_grab(void);
        void check_prog_name(char * prog_name, adios& my_adios);
        void check_comm_rank(char * comm_rank);
        void check_thread(std::string& thread);
};

}; // end namespace extractor
