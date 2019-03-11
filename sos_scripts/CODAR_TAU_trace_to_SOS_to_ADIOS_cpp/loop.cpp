/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "utils.hpp"
#include "sos_wrapper.hpp"
#include "adios_wrapper.hpp"
#include <nlohmann/json.hpp>

void main_loop(extractor::sos& my_sos, extractor::adios& my_adios) {
    PRINTSTACK
    bool done{false};
    int frame = 0;
    do {
        /* Check for complete new frame */
        if (my_sos.check_for_frame(frame)) {
            /* Write the new metadata */
            my_sos.write_metadata(frame, my_adios);
            /* Write the new events */
            my_sos.write_timer_data(frame, my_adios);
            frame++;
        } else { 
            done = true; 
        }
        /* If no more new frames, exit */
    } while (!done);
}
