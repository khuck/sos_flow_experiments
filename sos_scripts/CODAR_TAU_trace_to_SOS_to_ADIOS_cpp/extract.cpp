/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "sos.h"
#include "adios2.h"
#include <iostream>

#define DEBUG
#ifdef DEBUG
#define PRINTSTACK std::cout << __func__ << std::endl;
#else // DEBUG
#define PRINTSTACK 
#endif // DEBUG

namespace extractor {

/* Class containing SOS connection info */
class sos {
    private:
        bool connected;
    public:
        sos() : connected(false) {
            PRINTSTACK
        };
        ~sos() {
            PRINTSTACK
            if (connected) {
                this->disconnect();
            }
        };
        void disconnect() {
            PRINTSTACK
        };
};

/* Class containing ADIOS archive info */
class adios {
    private:
        bool opened;
    public:
        adios() : opened(false) {
            PRINTSTACK
        };
        ~adios() {
            PRINTSTACK
            if (opened) {
                this->close();
            }
        };
        void close() {
            PRINTSTACK
        };
};

}; // end namespace extractor

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