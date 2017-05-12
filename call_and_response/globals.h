#pragma once

#include <string>

#ifndef MY_EXTERN
#define MY_EXTERN extern
MY_EXTERN int _commrank;
MY_EXTERN int _commsize;
MY_EXTERN SOS_pub *_sos_pub;
MY_EXTERN SOS_runtime * _runtime;
MY_EXTERN int _daemon_rank;
MY_EXTERN bool _shutdown_daemon;
#endif //MY_EXTERN

void make_pub(void);
void initialize(int * argc, char *** argv);
void do_fork(std::string forkCommand);
void fork_exec_sosd_shutdown(void);
void send_shutdown_message(void);
void fork_exec_sosd(void);
void finalize(void);
void send_data(void);


