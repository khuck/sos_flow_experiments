#pragma once

#include "mpi.h"
#include <string>
#include <iostream>
#include <atomic>
#include "matrix.h"

/* Make sure to declare & initialize these in main.cpp! */
#ifndef MY_EXTERN
#define MY_EXTERN extern
MY_EXTERN int _commrank;
MY_EXTERN int _commsize;
MY_EXTERN int _daemon_rank;
MY_EXTERN bool _shutdown_daemon;
MY_EXTERN bool _balanced;
MY_EXTERN std::atomic<bool> _got_message;
#endif //MY_EXTERN

void make_pub(void);
void initialize(int * argc, char *** argv);
void do_fork(std::string forkCommand);
void send_shutdown_message(void);
void finalize(void);
void main_loop(int iterations, int write_iteration);
void setup_system_data(void);
void flush_it(void);
void do_broadcast(Matrix<double> &matrix);
double do_reduction(Matrix<double> &matrix);
void do_alltoall(Matrix<double> &matrix);
void do_adios(Matrix<double> &matrix);

#define zprint(...) { \
  if (_commrank == 0) { \
    printf(__VA_ARGS__); \
    fflush(stdout); \
  } \
}


