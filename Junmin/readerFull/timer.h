#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>


extern uint64_t _lastMeasuredMillis ;
//uint64_t _lastMeasuredMillis = 0;

uint64_t getCurrentTimeMillis();

#define logTimeMillis2(...) {}

#define  logTimeMillis(...) {    \
  uint64_t  ms = getCurrentTimeMillis();	\
  uint64_t  d = ms - _lastMeasuredMillis; \
  printf("   ELAPSED: %ld millis\t", d);\
  uint64_t  sec = d/1000; \
  if (sec > 0 ) \
    printf(" = %d sec \t", d/1000);\
  uint64_t  min = sec/60; \
  if (min > 0) \
    printf(" = %d min:%d sec ", min, (sec-min*60)); \
  printf(", "); \
  fprintf ( stdout, __VA_ARGS__);\
  fflush(stdout); \
  printf("\n"); \
  _lastMeasuredMillis = ms; \
  }

#endif
