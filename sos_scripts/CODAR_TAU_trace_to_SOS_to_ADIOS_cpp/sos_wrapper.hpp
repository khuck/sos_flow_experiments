/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#pragma once

#include "utils.hpp"
#include <nlohmann/json.hpp>
#include <string.h>
#include <unordered_map>
#include <vector>

using json = nlohmann::json;

namespace extractor {

/* Forward declaration, included in adios_wrapper.cpp */
extern class adios;

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
        int max_comm_rank;
        int max_threads;
        std::unordered_map<std::string, int> prog_names;
        std::unordered_map<std::string, int> value_names;
        std::unordered_map<std::string, int> metadata_keys;
        std::unordered_map<std::string, int> groups;
        std::unordered_map<std::string, int> timers;
        std::unordered_map<std::string, int> event_types;
        std::unordered_map<std::string, int> counters;
        static std::vector<std::string> split(const std::string& s, char delimiter);
        static bool replace(std::string& str, const std::string& from, const std::string& to);
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
            time_index(-1),
            max_comm_rank(0),
            max_threads(0)
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
        int check_prog_name(char * prog_name, adios& my_adios);
        int check_event_type(std::string& event_type, adios& my_adios);
        int check_timer(std::string& timer, adios& my_adios);
        int check_counter(std::string& timer, adios& my_adios);
        // do we need this?
        void check_comm_rank(int comm_rank);
        // do we need this?
        void check_thread(int thread);
        int get_prog_count(void) { return prog_names.size(); }
        int get_value_name_count(void) { return value_names.size(); }
        int get_timer_count(void) { return timers.size(); }
        int get_event_type_count(void) { return event_types.size(); }
        int get_counter_count(void) { return counters.size(); }
        int get_comm_rank_count(void) { return max_comm_rank+1; }
        int get_thread_count(void) { return max_threads+1; }
};

}; // end namespace extractor
