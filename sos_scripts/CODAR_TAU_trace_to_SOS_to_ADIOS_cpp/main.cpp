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
    TAU_SCOPED_TIMER_FUNC();

    std::string filename{"./config.json"};
    if (argc == 1) {
        // use the default
    } else if (argc == 2) {
        filename = std::string(argv[1]);
    } else if (argc > 3) {
        std::cout << "Usage: extract config.json" << std::endl;
    }
    std::cout << "Using config file: " << filename << std::endl;
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
