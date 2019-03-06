/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "utils.hpp"
#include "sos_wrapper.hpp"
#include "adios_wrapper.hpp"

int main (int argc, char * argv[]) {

    /* Connect to SOS */
    extractor::sos my_sos{};

    /* Open ADIOS archive */
    extractor::adios my_adios{};

    /* Loop until end of data */

    /* Close ADIOS archive */
    my_adios.close();

    /* Disconnect from SOS */
    my_sos.disconnect();

    return 0;
}