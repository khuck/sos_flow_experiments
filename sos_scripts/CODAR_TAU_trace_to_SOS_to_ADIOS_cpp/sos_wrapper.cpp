/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "utils.hpp"
#include "ssos.h"
#include "sosa.h"
#include <iostream>
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
        for (int r = 0 ; r < results.row_count ; r++) {
            if (atoi(results.data[r][pf_index]) < frame) {
                ready = false;
            }
        }
        SSOS_result_destroy(&results);
        if (!ready) { 
            if (use_timeout && (not_ready_loops > max_loops)) {
                quit = true;
                return false;
            }
            not_ready_loops++;
            std::cout << ".";
            fflush(stdout);
            usleep(config["sosd"]["usec_between_queries"].get<int>()); 
        }
    } while (!ready && !quit);
    std::cout << "Got frame " << frame << std::endl;
    return true;
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
}

}; // end namespace extractor
