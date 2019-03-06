/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#include "utils.hpp"
#include "adios2.h"

namespace extractor {

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
