/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "utils.hpp"
#include "sos_wrapper.hpp"
#include "adios_wrapper.hpp"
#include <nlohmann/json.hpp>

int main (int argc, char * argv[]) {

    std::string filename{"./config.json"};
    std::ifstream config_file(filename);
    json config;
    config_file >> config;
    //std::cout << std::setw(2) << config << std::endl;

    /* Connect to SOS */
    extractor::sos my_sos{config};

    /* Open ADIOS archive */
    extractor::adios my_adios{config};

    /* Loop until end of data */

    /* Close ADIOS archive */
    my_adios.close();

    /* Disconnect from SOS */
    my_sos.disconnect();

    return 0;
}