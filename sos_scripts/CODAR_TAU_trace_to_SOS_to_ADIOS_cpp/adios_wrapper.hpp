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
        adios2::Variable<size_t> timer_event_count;
        adios2::Variable<int> counter_count;
        adios2::Variable<size_t> counter_event_count;
        adios2::Variable<size_t> comm_count;
        adios2::Variable<unsigned long> event_timestamps;
        adios2::Variable<unsigned long> counter_values;
        adios2::Variable<unsigned long> comm_timestamps;
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
        void write_variables(sos& my_sos, size_t num_timer_values,
                size_t num_counter_values, size_t num_comm_values,
                std::vector<unsigned long>& timer_values_array,
                std::vector<unsigned long>& counter_values_array,
                std::vector<unsigned long>& comm_values_array);
};

}; // end namespace extractor
