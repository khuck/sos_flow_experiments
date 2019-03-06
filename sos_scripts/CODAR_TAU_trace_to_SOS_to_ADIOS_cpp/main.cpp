/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "utils.hpp"
#include "sos_wrapper.hpp"
#include "adios_wrapper.hpp"
#include <nlohmann/json.hpp>
#include <fstream>

// Defined in loop.cpp
void main_loop(extractor::sos& my_sos, extractor::adios& my_adios);

int main (int argc, char * argv[]) {
    PRINTSTACK

    std::string filename{"./config.json"};
    std::ifstream config_file(filename);
    json config;
    config_file >> config;

    /* Connect to SOS */
    extractor::sos my_sos{config};

    /* Open ADIOS archive */
    extractor::adios my_adios{config};

    /* Loop until end of data */
    main_loop(my_sos,  my_adios);

    /* Close ADIOS archive */
    my_adios.close();

    /* Disconnect from SOS */
    my_sos.disconnect();

    return 0;
}