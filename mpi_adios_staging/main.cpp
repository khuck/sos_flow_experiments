#define MY_EXTERN
#include "globals.h"
#include "pthread.h"
#include <thread>
#include <unistd.h>
#include "adios.h"

// instantiate globals.
int _commrank = 0;
int _commsize = 1;
int _daemon_rank = 0;
bool _shutdown_daemon = false;
bool _balanced = false;
bool _got_message = false;
const int _max_iterations = 10;
const int _write_iteration = 1;

inline unsigned int my_hardware_concurrency()
{
    unsigned int cores = std::thread::hardware_concurrency();
    return cores ? cores : sysconf(_SC_NPROCESSORS_ONLN);
}

void set_affinity() {
#if !defined(__APPLE__) && !defined(_MSC_VER)
  int s;
  cpu_set_t cpuset;
  pthread_t thread;
  thread = pthread_self();
  CPU_ZERO(&cpuset);
  unsigned int core = _commrank % my_hardware_concurrency();
  CPU_SET(core, &cpuset);

  //printf("rank %d binding to core %d.\n", _commrank, core);
  s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
  if (s != 0) {
    errno = s;
    perror("pthread_setaffinity_np");
    exit(EXIT_FAILURE);
  }
#endif // !defined(__APPLE__) && !defined(_MSC_VER)
}

int main (int argc, char * argv[]) {
  int iterations = _max_iterations;
  int write_iteration = _write_iteration;
  if (argc > 1) {
    iterations = atoi(argv[1]);
    write_iteration = atoi(argv[2]);
  }
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &_commrank);
  MPI_Comm_size(MPI_COMM_WORLD, &_commsize);
  set_affinity();
  adios_init ("matrix.xml", MPI_COMM_WORLD);
  printf("rank %d of %d, checking in.\n", _commrank, _commsize);
  MPI_Barrier(MPI_COMM_WORLD);
  main_loop(iterations, write_iteration);
  zprint("finalizing.\n");
  adios_finalize (_commrank);
  MPI_Finalize();
  return 0;
}
