/*
 * C++ program to extract TAU trace data from SOS listeners, and 
 * write it out to an ADIOS2 data file/stream.
 */

#pragma once

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// #define DEBUG
#ifdef DEBUG

namespace extractor {

class stackPrinter {
    private:
        static int depth;
    public:
        stackPrinter(const char * name) {
            for (int i = 0 ; i < depth ; i++) {
                std::cout << "-"; 
            } 
            std::cout << " " << name << std::endl;
            depth++;
        }
        ~stackPrinter() { depth--; }
};

inline int stackPrinter::depth{0};

}; // end namespace extractor

#define PRINTSTACK extractor::stackPrinter tmp{__func__};

#else // DEBUG
#define PRINTSTACK 
#endif // DEBUG

inline void _my_assert(const char* expression, const char* file, int line)
{
    fprintf(stderr, "Assertion '%s' failed, file '%s' line '%d'.",
        expression, file, line);
    abort();
}

#ifdef NDEBUG
#define MY_ASSERT(EXPRESSION) ((void)0)
#else
#define MY_ASSERT(EXPRESSION) ((EXPRESSION) ? (void)0 : \
    _my_assert(#EXPRESSION, __FILE__, __LINE__))
#endif
