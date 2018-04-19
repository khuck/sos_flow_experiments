#include "globals.h"
#include "simple_timer.h"
#include "matrix.h"
#include <algorithm>

#define MATRIX_SIZE 512
int NRA = MATRIX_SIZE; /* number of rows in matrix A */
int NCA = MATRIX_SIZE; /* number of columns in matrix A */
int NCB = MATRIX_SIZE; /* number of columns in matrix B */

void initialize(Matrix<double> &matrix) {
  int i,j;
#pragma omp parallel private(i,j) shared(matrix)
  {
    //set_num_threads();
    /*** Initialize matrices ***/
#pragma omp for nowait
    for (i=0; i<matrix.Rows(); i++) {
      for (j=0; j<matrix.Columns(); j++) {
        matrix[i][j]= i+j;
      }
    }
  }
}

// cols_a and rows_b are the same value
void compute(Matrix<double> &a, Matrix<double> &b, Matrix<double> &c) {
  simple_timer t("Compute");
  int i,j,k;
#pragma omp parallel private(i,j,k) shared(a,b,c)
  {
#pragma omp for nowait
    for (i=0; i<a.Rows(); i++) {
      for (k=0; k<a.Columns(); k++) {
        for(j=0; j<b.Columns(); j++) {
          c[i][j] += a[i][k] * b[k][j];
        }
      }
    }
  }   /*** End of parallel region ***/
}

double do_work(int i, bool do_adios_write) {
  //std::cout << _commrank << ": MATRIX SIZE: " << NRA << std::endl; fflush(stdout);

  Matrix<double> a(NRA, NCA);
  Matrix<double> b(NCA, NCB);
  Matrix<double> c(NRA, NCB);  

  /*** Spawn a parallel region explicitly scoping all variables ***/

  initialize(a);
  initialize(b);
  initialize(c);

  /* do a big broadcast */
  do_broadcast(a);
  do_broadcast(b);

  compute(a, b, c);

  /* do a reduction */
  double result = do_reduction(c);

  /* do an alltoall */
  do_alltoall(c);

  /* Do the ADIOS output */
  if (do_adios_write) {
    do_adios(c);
  }

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

void do_broadcast(Matrix<double> &matrix) {
  simple_timer t("Broadcast");
  int count = matrix.Rows() * matrix.Columns();
  MPI_Datatype datatype = MPI_DOUBLE;
  int root = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  void * buffer = matrix.Data();
  int rc = 0;
  mpi_check(MPI_Bcast( buffer, count, datatype, root, comm ));
}

double do_reduction(Matrix<double> &matrix) {
  simple_timer t("Reduction");
  int count = matrix.Rows() * matrix.Columns();
  MPI_Datatype datatype = MPI_DOUBLE;
  MPI_Op op = MPI_SUM;
  MPI_Comm comm = MPI_COMM_WORLD;
  double sum = 0.0;
  Matrix<double> sum_matrix(matrix.Rows(), matrix.Columns());
  const void * sendbuf = matrix.Data();
  void * recvbuf = sum_matrix.Data();
  mpi_check(MPI_Allreduce( sendbuf, recvbuf, count, datatype, op, comm ));
  // reduce the matrix, for fun.
  for (int i = 0 ; i < matrix.Rows() ; i++) {
    for (int j = 0 ; j < matrix.Columns() ; j++) {
      sum = sum + sum_matrix[i][j];
    }
  }
  return sum;
}

void do_alltoall(Matrix<double> &matrix) {
  simple_timer t("Alltoall");
  int count = matrix.Rows() * matrix.Columns();
  double * sendbuf = (double*)malloc((sizeof(double)) * count * _commsize);
  int index = _commrank * count;
  for (int i = 0 ; i < matrix.Rows() ; i++) {
    for (int j = 0 ; j < matrix.Columns() ; j++) {
      sendbuf[index++] = matrix[i][j];
    }
  }
  void * recvbuf = malloc((sizeof(double)) * count * _commsize);
  MPI_Datatype datatype = MPI_DOUBLE;
  MPI_Comm comm = MPI_COMM_WORLD;
  mpi_check(MPI_Alltoall( sendbuf, count, datatype, recvbuf, count, datatype, comm ));
}

void main_loop(int iterations, int write_iteration) {
  int i;
  double total = 0;
  simple_timer tm("Total Time");
  for (i = 0 ; i < iterations ; i++ ) {
    /* output status */
    if (_commrank == 0) {
      std::cout << "iteration " << i << std::endl; fflush(stdout);
    }
    /* wait for everyone to start at the same time */
    MPI_Barrier(MPI_COMM_WORLD);
    {
      simple_timer t("Iteration");
      /* do work */
      total += do_work(i, (((i+1)%write_iteration) == 0));
    }
    /* wait for everyone to finish at the same time */
    MPI_Barrier(MPI_COMM_WORLD);
  }
  /* make sure the final value is used */
  if (_commrank == 0) {
    std::cout << "Total: " << total << std::endl; fflush(stdout);
  }
}
