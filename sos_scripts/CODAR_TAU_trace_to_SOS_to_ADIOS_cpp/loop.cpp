/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "utils.hpp"
#include "sos_wrapper.hpp"
#include "adios_wrapper.hpp"
#include <nlohmann/json.hpp>
#include <unistd.h>

void main_loop(extractor::sos& my_sos, extractor::adios& my_adios) {
    PRINTSTACK
    bool done{false};
    int frame = 0;
    do {
        /* Check for complete new frame */
        my_sos.check_for_frame(frame);
        /* Query new frame of events */
        /* Write the new events */
        /* If no more new frames, exit */
        //done = true;
        frame++;
        sleep(1);
    } while (!done);
}
