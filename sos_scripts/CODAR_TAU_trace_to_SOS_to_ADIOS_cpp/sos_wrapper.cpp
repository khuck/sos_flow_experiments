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
#include <unistd.h>

namespace extractor {

void sos::connect() {
    PRINTSTACK
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
    PRINTSTACK
    if (connected) {
        SSOS_finalize();
    }
    connected = false;
};

bool sos::check_for_frame(int frame) {
    SSOS_query_results results;
    int max_frame_overall{0};
    bool ready = true;
    bool use_timeout = config["sosd"]["server_timeout"].get<bool>();
    int max_loops = config["sosd"]["exit_after_n_timeouts"].get<int>();
    int expected_pubs = config["sosd"]["aggregators"]["expected_pubs"].get<int>();
    int not_ready_loops = 0;
    bool quit = false;
    do {
        SSOS_request_pub_manifest(&results, &max_frame_overall, "", hostname.c_str(), portnumber);
        // std::cout << "Max frame: " << max_frame_overall << std::endl;
        // SOSA_results_output_to(stdout, reinterpret_cast<SOSA_results*>(&results), "", 3);
        /* find the "pub_frame" column */
        const char * pf = "pub_frame";
        int pf_index;
        for (int c = 0 ; c < results.col_count ; c++) {
            if (strcmp(results.col_names[c], pf) == 0) {
                pf_index = c;
                break;
            }
        }
        /* check all the pubs to make sure they are at the desired frame */
        if (results.row_count < expected_pubs) {
            ready = false;
            std::cout << expected_pubs << " pubs not seen yet." << std::endl;
        } else {
            for (int r = 0 ; r < results.row_count ; r++) {
                if (atoi(results.data[r][pf_index]) < frame) {
                    ready = false;
                    std::cout << "pubs not at frame " << frame << " yet." << std::endl;
                }
            }
        }
        SSOS_result_destroy(&results);
        if (!ready) { 
            if (use_timeout && (not_ready_loops > max_loops)) {
                quit = true;
                return false;
            }
            not_ready_loops++;
            usleep(config["sosd"]["usec_between_queries"].get<int>()); 
        }
    } while (!ready && !quit);
    std::cout << "Got frame " << frame << std::endl;
    static bool first_frame{true};
    if (first_frame) {
        build_column_map();
        first_frame = false;
    }
    return true;
}

void sos::build_column_map() {
    SSOS_query_results results;
    SSOS_cache_grab("", "", 0, 1, hostname.c_str(), portnumber);
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
    /*
    for (auto element : column_map) {
        std::cout << element.first << "::" << element.second << std::endl;
    }
    */
}

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

void sos::check_prog_name(char * prog_name, adios& my_adios) {
    if (prog_names.count(prog_name) == 0) {
        std::stringstream ss;
        ss << "program_name " << prog_names.size();
        prog_names[prog_name] = prog_names.size();
        my_adios.define_attribute(ss.str(), prog_name);
    }
}

void sos::check_comm_rank(char * comm_rank) {
    if (comm_ranks.count(comm_rank) == 0) {
        comm_ranks[comm_rank] = comm_ranks.size();
    }
}

void sos::check_thread(std::string& thread) {
    if (threads.count(thread) == 0) {
        threads[thread] = threads.size();
    }
}

void sos::write_metadata(int frame, adios& my_adios) {
    SSOS_query_results results;
    SSOS_cache_grab("", "TAU_Metadata", frame, 1, hostname.c_str(), portnumber);
    SSOS_result_claim(&results);
    //SOSA_results_output_to(stdout, reinterpret_cast<SOSA_results*>(&results), "", 0);
    int total_valid = results.row_count;
    // iterate over the rows
    for (int r = 0 ; r < results.row_count ; r++) {
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
        check_comm_rank(comm_rank);
        // tease apart the metadata name.
        auto tokens = split(value_name, ':');
        // thread = tokens[1]
        // metadata_key = tokens[2]
        check_thread(tokens[1]);
        std::stringstream attr_name;
        attr_name << "MetaData:" << prog_names[prog_name];
        attr_name << ":" << comm_rank << ":" << tokens[1] << ":" << tokens[2];
        my_adios.define_attribute(attr_name.str(), value);
    }
    SSOS_result_destroy(&results);
}

void sos::write_timer_data(int frame, adios& my_adios) {
    SSOS_query_results results;
    SSOS_cache_grab("", "TAU_EVENT", frame, 1, hostname.c_str(), portnumber);
    SSOS_result_claim(&results);
    total_valid = results.row_count;
    // SOSA_results_output_to(stdout, reinterpret_cast<SOSA_results*>(&results), "", 0);
    std::vector<long> timer_values_array{total_valid * 6};
    std::vector<long> counter_values_array{total_valid * 6};
    std::vector<long> comm_values_array{total_valid * 8};

    int timer_index = 0;
    int counter_index = 0;
    int comm_index = 0;
    /* create sorted tree of indexes into the results, ordered by
     * the timestamp, which is the value */
    std::vector< std::pair <long,int> > sorted;
    total_valid = results.row_count;
    for (int r = 0 ; r < results.row_count ; r++) {
        // can't use the timestamp, must sort by the pack time stamp.
        //sorted.push_back(std::make_pair(atol(results.data[r][value_index]),r));
        sorted.push_back(std::make_pair((atof(results.data[r][time_index]) * 1000000),r));
    }
    sort(sorted.begin(), sorted.end());
    for(auto iter = sorted.begin(); iter != sorted.end(); ++iter) {
        long timestamp = iter->first;
        int r = iter->second;
        std::cout << timestamp << " " << results.data[r][value_name_index] << std::endl;
        int row_frame = atoi(results.data[r][frame_index]);
        if (row_frame != frame) {
            total_valid = total_valid - 1;
            continue;
        }
        char * prog_name = results.data[r][prog_name_index];
        char * comm_rank = results.data[r][comm_rank_index];
        char * value_name = results.data[r][value_name_index];
        char * value = results.data[r][value_index];
        check_prog_name(prog_name, my_adios);
        check_comm_rank(comm_rank);
    }

    SSOS_result_destroy(&results);
}

void sos::sql_query() {
    SSOS_query_results results;
    SSOS_query_exec("select * from tblpubs;", hostname.c_str(), portnumber);
    SSOS_result_claim(&results);
    SOSA_results_output_to(stdout, reinterpret_cast<SOSA_results*>(&results), "", 0);
    SSOS_result_destroy(&results);
}

void sos::cache_grab() {
    SSOS_query_results results;
    SSOS_cache_grab("", "", -1, 10, hostname.c_str(), portnumber);
    SSOS_result_claim(&results);
    SOSA_results_output_to(stdout, reinterpret_cast<SOSA_results*>(&results), "", 0);
    SSOS_result_destroy(&results);
}

}; // end namespace extractor
