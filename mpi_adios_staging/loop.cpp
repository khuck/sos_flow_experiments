#include "globals.h"
#include "simple_timer.h"
#include <algorithm>

#define MATRIX_SIZE 1024
const int max_iterations = 10;
const double increment_divisor = 0.01; /* add/subtract 1% each time */
int NRA = MATRIX_SIZE; /* number of rows in matrix A */
int NCA = MATRIX_SIZE; /* number of columns in matrix A */
int NCB = MATRIX_SIZE; /* number of columns in matrix B */

double** allocateMatrix(int rows, int cols) {
  int i;
  double **matrix = (double**)malloc((sizeof(double*)) * rows);
  for (i=0; i<rows; i++) {
    matrix[i] = (double*)malloc((sizeof(double)) * cols);
  }
  return matrix;
}

void freeMatrix(double** matrix, int rows, int cols) {
  int i;
  for (i=0; i<rows; i++) {
    free(matrix[i]);
  }
  free(matrix);
}

void initialize(double **matrix, int rows, int cols) {
  int i,j;
#pragma omp parallel private(i,j) shared(matrix)
  {
    //set_num_threads();
    /*** Initialize matrices ***/
#pragma omp for nowait
    for (i=0; i<rows; i++) {
      for (j=0; j<cols; j++) {
        matrix[i][j]= i+j;
      }
    }
  }
}

// cols_a and rows_b are the same value
void compute(double **a, double **b, double **c, int rows_a, int cols_a, int cols_b) {
  int i,j,k;
#pragma omp parallel private(i,j,k) shared(a,b,c)
  {
#pragma omp for nowait
    for (i=0; i<rows_a; i++) {
      for (k=0; k<cols_a; k++) {
        for(j=0; j<cols_b; j++) {
          c[i][j] += a[i][k] * b[k][j];
        }
      }
    }
  }   /*** End of parallel region ***/
}

#define increment int((MATRIX_SIZE)*increment_divisor)

double do_work(int i) {
  double **a,    /* matrix A to be multiplied */
  **b,           /* matrix B to be multiplied */
  **c;           /* result matrix C */
  NRA = MATRIX_SIZE;
  NCA = MATRIX_SIZE;
  NCB = MATRIX_SIZE;
  //std::cout << _commrank << ": MATRIX SIZE: " << NRA << std::endl; fflush(stdout);

  a = allocateMatrix(NRA, NCA);
  b = allocateMatrix(NCA, NCB);
  c = allocateMatrix(NRA, NCB);  

  /*** Spawn a parallel region explicitly scoping all variables ***/

  initialize(a, NRA, NCA);
  initialize(b, NCA, NCB);
  initialize(c, NRA, NCB);

  /* do a big broadcast */
  do_broadcast(a, NRA, NCA);
  do_broadcast(b, NRA, NCA);

  compute(a, b, c, NRA, NCA, NCB);

  /* do a reduction */
  double result = do_reduction(c, NRA, NCB);

  /* do an alltoall */
  do_alltoall(c, NRA, NCA);

  freeMatrix(a, NRA, NCA);
  freeMatrix(b, NCA, NCB);
  freeMatrix(c, NCA, NCB);

  return result;
}

bool mpi_check(int rc) {
  switch (rc) {
    case MPI_SUCCESS:
		return true;
    case MPI_ERR_COMM:
        std::cerr << "Invalid communicator.\n";
        std::cerr << "A common error is to use a null communicator in a call \n";
        std::cerr << "(not even allowed in MPI_Comm_rank ).";
        std::cerr << std::endl;
        abort();
    case MPI_ERR_COUNT:
        std::cerr << "Invalid  count argument.\n";
        std::cerr << "Count arguments must be non-negative; a count of zero is often valid.";
        std::cerr << std::endl;
        abort();
    case MPI_ERR_TYPE:
        std::cerr << "Invalid datatype argument.\n";
        std::cerr << "Additionally, this error can occur  if  an  uncommitted MPI_Datatype\n";
        std::cerr << "(see MPI_Type_commit ) is used in a communication call.";
        std::cerr << std::endl;
        abort();
    case MPI_ERR_BUFFER:
        std::cerr << "Invalid buffer pointer.\n";
        std::cerr << "Usually a null buffer where one is not valid.";
        std::cerr << std::endl;
        abort();
    case MPI_ERR_ROOT:
        std::cerr << "Invalid  root.\n";
        std::cerr << "The  root must be specified as a rank in the communicator.\n";
        std::cerr << "Ranks must be between zero and the size of the communicator minus one.";
        std::cerr << std::endl;
        abort();
    default:
        return false;
  }
}

void do_broadcast(double ** matrix, int rows, int cols) {
  int count = cols;
  MPI_Datatype datatype = MPI_DOUBLE;
  int root = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  void * buffer = NULL;
  int rc = 0;
  for (int i = 0 ; i < rows ; i++) {
    buffer = &(matrix[i][0]);
    mpi_check(MPI_Bcast( buffer, count, datatype, root, comm ));
  }
}

double do_reduction(double ** matrix, int rows, int cols) {
  const void * sendbuf;
  int count = cols;
  MPI_Datatype datatype = MPI_DOUBLE;
  MPI_Op op = MPI_SUM;
  MPI_Comm comm = MPI_COMM_WORLD;
  double sum = 0.0;
  double * recvbuf = (double*)malloc((sizeof(double)) * cols);
  for (int i = 0 ; i < rows ; i++) {
    sendbuf = &(matrix[i][0]);
    mpi_check(MPI_Allreduce( sendbuf, (void*)recvbuf, count, datatype, op, comm ));
    for (int j = 0 ; j < cols ; j++) {
      sum = sum + recvbuf[j];
    }
  }
  free(recvbuf);
  return sum;
}

void do_alltoall(double ** matrix, int rows, int cols) {
}

void main_loop(void) {
  int i;
  double total = 0;
  simple_timer tm("Total Time");
  for (i = 0 ; i < max_iterations ; i++ ) {
    /* output status */
    if (_commrank == 0) {
      std::cout << "iteration " << i << std::endl; fflush(stdout);
    }
    /* wait for everyone to start at the same time */
    MPI_Barrier(MPI_COMM_WORLD);
    {
      simple_timer t("Iteration");
      /* do work */
      total += do_work(i);
    }
    /* wait for everyone to finish at the same time */
    MPI_Barrier(MPI_COMM_WORLD);
  }
  /* make sure the final value is used */
  if (_commrank == 0) {
    std::cout << "Total: " << total << std::endl; fflush(stdout);
  }
}
