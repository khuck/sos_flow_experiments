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

/* Forward declaration, included in adios_wrapper.cpp */
extern class sos;

/* Class containing ADIOS archive info */
class adios {
    private:
        bool opened;
        json config;
        adios2::ADIOS ad;
        adios2::IO bpIO;
        adios2::Engine bpWriter;
        adios2::Variable<int> program_count;
        adios2::Variable<int> comm_rank_count;
        adios2::Variable<int> thread_count;
        adios2::Variable<int> event_type_count;
        adios2::Variable<int> timer_count;
        adios2::Variable<int> timer_event_count;
        adios2::Variable<int> counter_count;
        adios2::Variable<int> counter_event_count;
        adios2::Variable<int> comm_count;
        adios2::Variable<long> event_timestamps;
        adios2::Variable<long> counter_values;
        adios2::Variable<long> comm_timestamps;
    public:
        adios(json& _config) : 
            opened(false),
            config(_config)
        {
            PRINTSTACK
            initialize();
            open();
            define_variables();
        };
        ~adios() {
            PRINTSTACK
            close();
        };
        void initialize();
        void define_variables();
        void open();
        void close();
        void define_attribute(std::string name, std::string value);
        void write_variables(sos& my_sos, int num_timer_values,
                int num_counter_values, int num_comm_values,
                std::vector<long>& timer_values_array,
                std::vector<long>& counter_values_array,
                std::vector<long>& comm_values_array);
};

}; // end namespace extractor
