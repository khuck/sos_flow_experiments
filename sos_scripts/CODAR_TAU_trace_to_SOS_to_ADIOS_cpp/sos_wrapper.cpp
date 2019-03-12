/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "utils.hpp"
#include "ssos.h"
#include "sosa.h"
#include <iostream>
#include <vector>
#include <string.h>
#include "sos_wrapper.hpp"
#include "adios_wrapper.hpp"
#include <unistd.h>
#include "stdlib.h"

namespace extractor {

    /* Connect to the SOS daemon using the hostname and port specified
     */
    void sos::connect() {
        SSOS_init("extractor");
        connected = true;
        if (config["sosd"]["SOS_CMD_PORT"] != nullptr) {
            portnumber = config["sosd"]["SOS_CMD_PORT"].get<int>();
        }
        if (config["sosd"]["hostname"] != nullptr) {
            hostname = config["sosd"]["hostname"].get<std::string>();
        }
        test_connection();
    }

    void sos::test_connection() {
        int online{0};
        SSOS_is_online(&online);
        if (online == 0) {
            std::cerr << "SOS not online" << std::endl;
            exit(9);
        }
    }

    void sos::disconnect() {
        SSOS_result_destroy(&results);
        if (connected) {
            SSOS_finalize();
        }
        connected = false;
    };

    /* Get the SOS manifest, and see how many pubs have reached/passed
     * the frame of interest.
     */
    bool sos::check_for_frame(int frame) {
        TAU_SCOPED_TIMER_FUNC();
        static SSOS_query_results * manifest_results;
        static bool first{true};
        int max_frame_overall{0};
        bool ready = true;
        bool use_timeout = config["sosd"]["server_timeout"].get<bool>();
        int max_loops = config["sosd"]["exit_after_n_timeouts"].get<int>();
        int expected_pubs = config["sosd"]["aggregators"]["expected_pubs"].get<int>();
        int not_ready_loops = 0;
        bool quit = false;
        do {
            ready = true;
            {
                /* The first time we get the manifest, it has to be initialized.
                 * After that, we can just refresh the values to avoid memory
                 * free/alloc cycles
                 */
                TAU_SCOPED_TIMER("SSOS_request_pub_manifest");
                if (UNLIKELY(first)) {
                    SSOS_request_pub_manifest(&manifest_results, &max_frame_overall, "", hostname.c_str(), portnumber);
                    first = false;
                } else {
                    SSOS_refresh_pub_manifest(manifest_results, &max_frame_overall, "", hostname.c_str(), portnumber);
                }
            }
            /* find the "pub_frame" column */
            const char * pf = "pub_frame";
            static int pf_index = -1;
            if (UNLIKELY(pf_index == -1)) {
                for (int c = 0 ; c < manifest_results->col_count ; c++) {
                    if (strcmp(manifest_results->col_names[c], pf) == 0) {
                        pf_index = c;
                        break;
                    }
                }
            }
            /* check all the pubs to make sure they are at the desired frame */
            if (manifest_results->row_count < expected_pubs) {
                ready = false;
                std::cout << expected_pubs << " pubs not seen yet." << std::endl;
            } else {
                for (int r = 0 ; r < manifest_results->row_count ; r++) {
                    if (atoi(manifest_results->data[r][pf_index]) < frame) {
                        ready = false;
                        std::cout << "pubs not at frame " << frame << " yet." << std::endl;
                    }
                }
            }
            /* If we are still waiting, count how many timeouts we have,
             * and don't exceed the number specified in the config file
             */
            if (!ready) { 
                not_ready_loops++;
                if (use_timeout && (not_ready_loops >= max_loops)) {
                    quit = true;
                    std::cout << "Exiting main loop." << std::endl;
                    /* No more work, exit. */
                    return false;
                }
                usleep(config["sosd"]["usec_between_queries"].get<int>()); 
            }
        } while (!ready && !quit);
        /* All pubs are at the frame, we are good to go! */
        std::cout << "Got frame " << frame << std::endl;
        /* First, get the column names for the cache grabs */
        static bool first_frame{true};
        if (UNLIKELY(first_frame)) {
            build_column_map();
            first_frame = false;
        }
        return true;
    }

    /* Do a cache grab, and determine what order the columns are in */
    void sos::build_column_map() {
        TAU_SCOPED_TIMER_FUNC();
        {
            TAU_SCOPED_TIMER("SSOS_cache_grab");
            SSOS_cache_grab("", "", 0, 1, hostname.c_str(), portnumber);
        }
        // this is the first claim, so it will initialize the results
        SSOS_result_claim(&results);
        for (int c = 0 ; c < results.col_count ; c++) {
            column_map[results.col_names[c]] = c;
        }
        prog_name_index = column_map["prog_name"];
        comm_rank_index = column_map["comm_rank"];
        value_name_index = column_map["val_name"];
        value_index = column_map["val"];
        frame_index = column_map["frame"];
        time_index = column_map["time_pack"];
    }

    /* Helper function to tokenize strings into vectors of strings */
    std::vector<std::string> sos::split(const std::string& s, char delimiter)
    {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream tokenStream(s);
        while (std::getline(tokenStream, token, delimiter))
        {
            tokens.push_back(token);
        }
        return tokens;
    }

    /* Helper function to do substring replacement */
    bool sos::replace(std::string& str, const std::string& from, const std::string& to) {
        size_t start_pos = str.find(from);
        if(start_pos == std::string::npos)
            return false;
        str.replace(start_pos, from.length(), to);
        return true;
    }

    /* Keep a map of program names to indexes */
    int sos::check_prog_name(char * prog_name, adios& my_adios) {
        if (prog_names.count(prog_name) == 0) {
            std::stringstream ss;
            ss << "program_name " << prog_names.size();
            prog_names[prog_name] = prog_names.size();
            my_adios.define_attribute(ss.str(), prog_name);
        }
        return prog_names[prog_name];
    }

    /* Keep a map of event types to indexes */
    int sos::check_event_type(std::string& event_type, adios& my_adios) {
        if (event_types.count(event_type) == 0) {
            std::stringstream ss;
            ss << "event_type " << event_types.size();
            event_types[event_type] = event_types.size();
            my_adios.define_attribute(ss.str(), event_type);
        }
        return event_types[event_type];
    }

    /* Keep a map of timers to indexes */
    int sos::check_timer(std::string& timer, adios& my_adios) {
        if (timers.count(timer) == 0) {
            std::stringstream ss;
            ss << "timer " << timers.size();
            timers[timer] = timers.size();
            my_adios.define_attribute(ss.str(), timer);
        }
        return timers[timer];
    }

    /* Keep a map of counters to indexes */
    int sos::check_counter(std::string& counter, adios& my_adios) {
        if (counters.count(counter) == 0) {
            std::stringstream ss;
            ss << "counter " << counters.size();
            counters[counter] = counters.size();
            my_adios.define_attribute(ss.str(), counter);
        }
        return counters[counter];
    }

    /* Keep a map of max number of MPI communicators */
    void sos::check_comm_rank(int comm_rank) {
        if (comm_rank > max_comm_rank) {
            max_comm_rank = comm_rank;
        }
    }

    /* Keep a map of max number of threads per process */
    void sos::check_thread(int thread) {
        if (thread > max_threads) {
            max_threads = thread;
        }
    }

    /* Write the TAU metadata to the trace as ADIOS attributes */
    void sos::write_metadata(int frame, adios& my_adios) {
        TAU_SCOPED_TIMER_FUNC();
        {
            TAU_SCOPED_TIMER("SSOS_cache_grab");
            SSOS_cache_grab("", "TAU_Metadata", frame, 1, hostname.c_str(), portnumber);
        }
        {
            TAU_SCOPED_TIMER("SSOS_result_claim");
            SSOS_result_claim_initialized(&results, 0);
        }
        int total_valid = results.row_count;
        /* iterate over the rows, order doesn't matter */
        for (int r = 0 ; r < results.row_count ; r++) {
            TAU_SCOPED_TIMER("sos::write_metadata for loop");
            int this_frame = atoi(results.data[r][frame_index]);
            if (this_frame != frame) {
                total_valid = total_valid - 1;
                continue;
            }
            char * prog_name = results.data[r][prog_name_index];
            char * comm_rank = results.data[r][comm_rank_index];
            char * value_name = results.data[r][value_name_index];
            char * value = results.data[r][value_index];
            if (strcmp(value, "") == 0) {
                continue;
            }
            check_prog_name(prog_name, my_adios);
            check_comm_rank(atoi(comm_rank));
            // tease apart the metadata name.
            auto tokens = split(value_name, ':');
            // thread = tokens[1]
            // metadata_key = tokens[2]
            check_thread(stoi(tokens[1]));
            std::stringstream attr_name;
            attr_name << "MetaData:" << prog_names[prog_name];
            attr_name << ":" << comm_rank << ":" << tokens[1] << ":" << tokens[2];
            my_adios.define_attribute(attr_name.str(), value);
        }
    }

    /* Get all the events in the last frame of SOS data
     * and write them to the trace as ADIOS variables */
    void sos::write_timer_data(int frame, adios& my_adios) {
        TAU_SCOPED_TIMER_FUNC();
        {
            TAU_SCOPED_TIMER("SSOS_cache_grab");
            SSOS_cache_grab("", "TAU_EVENT", frame, 1, hostname.c_str(), portnumber);
        }
        {
            TAU_SCOPED_TIMER("SSOS_result_claim");
            SSOS_result_claim_initialized(&results, 0);
        }
        TAU_START("overheads")
        total_valid = results.row_count;
        std::vector<unsigned long> timer_values_array(total_valid * 6, 0);
        std::vector<unsigned long> counter_values_array(total_valid * 6, 0);
        std::vector<unsigned long> comm_values_array(total_valid * 8, 0);

        int timer_value_index = 0;
        int counter_value_index = 0;
        int comm_value_index = 0;
        /* create sorted tree of indexes into the results, ordered by
         * the timestamp, which is the pack time (not the true timestamp).
         * This way, we ensure everything is in order for downstream processing */
        std::vector< std::pair <unsigned long,int> > sorted;
        total_valid = results.row_count;
        if (results.row_count == 0) {
            std::cerr << "Warning - no events in frame " << frame << std::endl;
        }
        for (int r = 0 ; r < results.row_count ; r++) {
            /* can't use the timestamp, must sort by the pack time stamp! */
            //sorted.push_back(std::make_pair(atol(results.data[r][value_index]),r));
            sorted.push_back(std::make_pair((atof(results.data[r][time_index]) * 1000000),r));
        }
        sort(sorted.begin(), sorted.end());
        static const std::string event_prefix{"TAU_EVENT_"};
        static const std::string empty{""};
        TAU_STOP("overheads")
        /* Iterate over the sorted timestamps, and process each 
         * associated row */
        for(auto iter = sorted.begin(); iter != sorted.end(); ++iter) {
            TAU_SCOPED_TIMER("sos::write_timer_data for loop");
            unsigned long timestamp = iter->first;
            int r = iter->second;
            int row_frame = atoi(results.data[r][frame_index]);
            if (row_frame != frame) {
                total_valid = total_valid - 1;
                continue;
            }
            char * prog_name = results.data[r][prog_name_index];
            int comm_rank = atoi(results.data[r][comm_rank_index]);
            char * value_name = results.data[r][value_name_index];
            char * value = results.data[r][value_index];
            int prog_index = check_prog_name(prog_name, my_adios);
            check_comm_rank(comm_rank);
            /* Handle timer start/stop events */
            if ((strstr(value_name, "TAU_EVENT_ENTRY") != NULL) ||
                    (strstr(value_name, "TAU_EVENT_EXIT") != NULL)) {
                // tease apart the event name.
                auto tokens = split(value_name, ':');
                replace(tokens[0], event_prefix, empty);
                int event_index = check_event_type(tokens[0], my_adios);
                int thread = stoi(tokens[1]);
                check_thread(thread);
                int timer_index = check_timer(tokens[2], my_adios);
                timer_values_array[timer_value_index++] = (unsigned long)(prog_index);
                timer_values_array[timer_value_index++] = (unsigned long)(comm_rank);
                timer_values_array[timer_value_index++] = (unsigned long)(thread);
                timer_values_array[timer_value_index++] = (unsigned long)(event_index);
                timer_values_array[timer_value_index++] = (unsigned long)(timer_index);
                timer_values_array[timer_value_index++] = timestamp;
            /* Handle counter events */
            } else if (strstr(value_name, "TAU_EVENT_COUNTER") != NULL) {
                // tease apart the counter name.
                auto tokens = split(value_name, ':');
                int thread = stoi(tokens[1]);
                check_thread(thread);
                int counter_index = check_counter(tokens[2], my_adios);
                counter_values_array[counter_value_index++] = (unsigned long)(prog_index);
                counter_values_array[counter_value_index++] = (unsigned long)(comm_rank);
                counter_values_array[counter_value_index++] = (unsigned long)(thread);
                counter_values_array[counter_value_index++] = (unsigned long)(counter_index);
                counter_values_array[counter_value_index++] = (unsigned long)(value);
                counter_values_array[counter_value_index++] = timestamp;
            /* Handle communication events */
            } else if ((strstr(value_name, "TAU_EVENT_SEND") != NULL) ||
                    (strstr(value_name, "TAU_EVENT_RECV") != NULL)) {
                // tease apart the counter name.
                auto tokens = split(value_name, ':');
                replace(tokens[0], event_prefix, empty);
                int event_index = check_event_type(tokens[0], my_adios);
                int thread = stoi(tokens[1]);
                check_thread(thread);
                counter_values_array[comm_value_index++] = (unsigned long)(prog_index);
                counter_values_array[comm_value_index++] = (unsigned long)(comm_rank);
                counter_values_array[comm_value_index++] = (unsigned long)(thread);
                counter_values_array[comm_value_index++] = (unsigned long)(event_index);
                counter_values_array[comm_value_index++] = atol(tokens[2].c_str());
                counter_values_array[comm_value_index++] = atol(tokens[3].c_str());
                counter_values_array[comm_value_index++] = atol(tokens[4].c_str());
                counter_values_array[comm_value_index++] = timestamp;
            } else {
                std::cerr << "Error: unknown event: " 
                    << prog_name << ", " 
                    << comm_rank << ", " 
                    << value_name << std::endl;
            }
        }
        /* Write to the ADIOS archive */
        my_adios.write_variables(*this, 
               timer_value_index/6, counter_value_index/6, comm_value_index/8,
               timer_values_array, counter_values_array, comm_values_array);
    }

    /* Reference */
    void sos::sql_query() {
        SSOS_query_exec("select * from tblpubs;", hostname.c_str(), portnumber);
        SSOS_result_claim(&results);
        SOSA_results_output_to(stdout, reinterpret_cast<SOSA_results*>(&results), "", 0);
    }

    /* Reference */
    void sos::cache_grab() {
        SSOS_cache_grab("", "", -1, 10, hostname.c_str(), portnumber);
        SSOS_result_claim(&results);
        SOSA_results_output_to(stdout, reinterpret_cast<SOSA_results*>(&results), "", 0);
    }

}; // end namespace extractor
