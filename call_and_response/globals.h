#pragma once

#include "mpi.h"
#include "sos.h"
#include <string>
#include <iostream>

/* Make sure to declare & initialize these in main.cpp! */
#ifndef MY_EXTERN
#define MY_EXTERN extern
MY_EXTERN int _commrank;
MY_EXTERN int _commsize;
MY_EXTERN SOS_pub *_sos_pub;
MY_EXTERN SOS_runtime * _runtime;
MY_EXTERN int _daemon_rank;
MY_EXTERN bool _shutdown_daemon;
MY_EXTERN bool _balanced;
#endif //MY_EXTERN

void make_pub(void);
void initialize(int * argc, char *** argv);
void do_fork(std::string forkCommand);
void fork_exec_sosd_shutdown(void);
void send_shutdown_message(void);
void fork_exec_sosd(void);
void finalize(void);
void sample_value(std::string, double);
void send_sos_system_data(void);
void main_loop(void);
void setup_system_data(void);

