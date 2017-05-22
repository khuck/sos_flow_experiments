#pragma once

#include <chrono>
#include <iostream>
#include <sstream>

class simple_timer {
    const double nanoseconds  = 1.0e9;
    const double microseconds = 1.0e6;
    const double milliseconds = 1.0e3;
    const double seconds    = 1.0;
  public:
    std::chrono::high_resolution_clock::time_point start;
    std::string _name;
    simple_timer() : start(std::chrono::high_resolution_clock::now()), _name("timer") {}
    simple_timer(std::string name) : start(std::chrono::high_resolution_clock::now()), _name(name) {}
    ~simple_timer() {
      std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - start);
      sample_value(_name, time_span.count() * seconds);
      /*
      std::stringstream ss;
      ss << _commrank << ":" << _name << ": " << time_span.count() * seconds << std::endl;
      std::cout << ss.str(); fflush(stdout);
      */
    }
};
