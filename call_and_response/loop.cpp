#include "globals.h"
#include "simple_timer.h"
#include <algorithm>

#define MATRIX_SIZE 512
const int max_iterations = 50;
const double increment_divisor = 1.0/max_iterations;
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
  double **a,           /* matrix A to be multiplied */
  **b,           /* matrix B to be multiplied */
  **c;           /* result matrix C */

  // upper ranks get an increase in work each iteration
  if ((i > 0) && (!_balanced)) {
    if (_commrank >= (_commsize / 2)) {
      NRA = std::min(int(MATRIX_SIZE*1.25),(NRA + increment));
      NCA = std::min(int(MATRIX_SIZE*1.25),(NCA + increment));
      NCB = std::min(int(MATRIX_SIZE*1.25),(NCB + increment));
    } else {
      NRA = std::max(int(MATRIX_SIZE*0.75),(NRA - increment));
      NCA = std::max(int(MATRIX_SIZE*0.75),(NCA - increment));
      NCB = std::max(int(MATRIX_SIZE*0.75),(NCB - increment));
    }
  } else {
    NRA = MATRIX_SIZE;
    NCA = MATRIX_SIZE;
    NCB = MATRIX_SIZE;
  }
  //std::cout << _commrank << ": MATRIX SIZE: " << NRA << std::endl;

  a = allocateMatrix(NRA, NCA);
  b = allocateMatrix(NCA, NCB);
  c = allocateMatrix(NRA, NCB);  

/*** Spawn a parallel region explicitly scoping all variables ***/

  initialize(a, NRA, NCA);
  initialize(b, NCA, NCB);
  initialize(c, NRA, NCB);

  compute(a, b, c, NRA, NCA, NCB);

  double result = c[0][1];

  freeMatrix(a, NRA, NCA);
  freeMatrix(b, NCA, NCB);
  freeMatrix(c, NCA, NCB);

  return result;
}

/* for now, just balance at 50% of the iterations complete.
 * Eventually, we want to get an analysis result from SOS. */
void check_for_balance(int i) {
    /* make the app balanced */
    if (i >= max_iterations/2) _balanced=true;
}

void main_loop(void) {
  int i;
  double total = 0;
  setup_system_data();
  for (i = 0 ; i < max_iterations ; i++ ) {
    /* Ask SOS for update */
    check_for_balance(i);
    /* make a timer */
    simple_timer t("Iteration");
    /* output status */
    if (_commrank == 0) {
      std::cout << "iteration " << i << std::endl; fflush(stdout);
    }
    /* wait for everyone to start at the same time */
    MPI_Barrier(MPI_COMM_WORLD);
    /* do work */
    total += do_work(i);
    /* output some system data */
    send_sos_system_data();
  }
  /* wait for everyone to finish at the same time */
  MPI_Barrier(MPI_COMM_WORLD);
  /* make sure the final value is used */
  if (_commrank == 0) {
    std::cout << "Total: " << total << std::endl;
  }
}