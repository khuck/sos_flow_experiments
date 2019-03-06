/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#pragma once

#include <iostream>

#define DEBUG
#ifdef DEBUG

namespace extractor {

class stackPrinter {
    private:
        static int depth;
    public:
        stackPrinter(const char * name) {
            for (int i = 0 ; i < depth ; i++) {
                std::cout << "--"; 
            } 
            std::cout << " " << name << std::endl;
            depth++;
        }
        ~stackPrinter() { depth--; }
};
int stackPrinter::depth{1};

}; // end namespace extractor

#define PRINTSTACK extractor::stackPrinter tmp{__func__};

#else // DEBUG
#define PRINTSTACK 
#endif // DEBUG

