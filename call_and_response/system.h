#pragma once

#include <vector>

class CPUStat {
public:
  char name[32];
  long long user;
  long long nice;
  long long system;
  long long idle;
  long long iowait;
  long long irq;
  long long softirq;
  long long steal;
  long long guest;
};

typedef std::vector<CPUStat*> CPUs;

class ProcData {
public:
  CPUs cpus;
  long long ctxt;
  long long btime;
  long processes;
  long procs_running;
  long procs_blocked;
  long power;
  long energy;
  long freshness;
  long generation;
  long power_cap;
  long startup;
  long version;
  long long package0;
  long long dram;
  ~ProcData(void);
  ProcData* diff(const ProcData& rhs);
  void sample_values();
};


