#include "mpi.h"
#include "sos.h"

#define MY_EXTERN
#include "globals.h"

// instantiate globals.
int _commrank = 0;
int _commsize = 1;
SOS_pub * _sos_pub = NULL;
SOS_runtime * _runtime = NULL;
int _daemon_rank = 0;
bool _shutdown_daemon = false;

int main (int argc, char * argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &_commrank);
  MPI_Comm_size(MPI_COMM_WORLD, &_commsize);
  printf("rank %d of %d, checking in.\n", _commrank, _commsize);
  initialize(&argc, &argv);
  send_data();
  finalize();
  MPI_Finalize();
  return 0;
}
