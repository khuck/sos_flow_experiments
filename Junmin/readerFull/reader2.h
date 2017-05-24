#ifndef READER_H
#define READER_H
/**
 *
 */

#include "mpi.h"
#include "adios.h"
#include "adios_read.h"

//#include "misc.h"
//#include "utils.h"
//#include "test_common.h"
//#include "cfg.h"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#include "timer.h"
// for printing the values of the variable
#define STR_BUFFER_SIZE 3000
#define NPOWER 4
#define NCOL 6

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


extern const char ATTR_NAMES[6] [20];

//typedef long double long_double_t ;
typedef  double long_double_t ;


struct _stats {
  long_double_t _max;
  long_double_t _min;
};

typedef struct _stats ARRAY_STAT;


struct _col {
  long_double_t* _colData;
  ARRAY_STAT _colStat;
};

typedef struct _col COL_DATA;


// 
// histogram
//
struct _Hist2dSpec {
  int _c1; // column id 0-6: x/Px/y/Py/zPz
  int _c2; // column id
  ARRAY_STAT _statC1; // all timesteps
  ARRAY_STAT _statC2; // all timesteps
  uint64_t _maxHits; // all timesteps
  char  _name1[20];
  char  _name2[20];
};

typedef struct _Hist2dSpec HISTOGRAM_2D_SPEC;

struct _var {
  long_double_t* _data;
  uint64_t _size;
  long_double_t _max;
  long_double_t _min;
};

typedef struct _var DATA;


struct _diag_stat {
  long_double_t _attr[NCOL][NPOWER];
  
  long_double_t _xpx;
  long_double_t _ypy;
  long_double_t _zpz;

  long_double_t _gam[2]; // gamma

  long_double_t _lrmax;

  long_double_t _maxByCol[NCOL];

  int _step;
  int _bunchID;
  uint64_t _bunchSize;
};

typedef struct _diag_stat DIAGNOSE;

// general
int getPos(const char* colName);


//diagnose
void clear(DIAGNOSE*);
void prepareDiagnose(COL_DATA* data, uint64_t nrows, int ncols, int rank, DIAGNOSE* d);
void finishDiagnose(DIAGNOSE* diagOfAllBunchs, int numBunch);

//
// 2D histogram related
//
void init2DHistInput(HISTOGRAM_2D_SPEC* input, int firstCol, int secondCol, const char* n1, const char* n2);
void intervine(int* h1, uint64_t h1size, int* h2, uint64_t h2size, int bin, uint64_t result[bin][bin]);
void binMap(DATA* v, int* result, int bin);
void calculate2DHistogram(DATA* v1, DATA* v2, int bin, uint64_t result[bin][bin], int rank);
long_double_t* getCol(long_double_t* data, uint64_t nrows, uint64_t ncols, uint64_t whichCol, ARRAY_STAT* stat);
void getAllCol(long_double_t* data, uint64_t nrows, uint64_t numCols, COL_DATA* colData);
//void process2D(long_double_t* data, uint64_t nrows, uint64_t ncols, int rank, int nbin, HISTOGRAM_2D_SPEC* input, int steps);
void process2D(COL_DATA* data, uint64_t nrows, uint64_t ncols, int rank, int nbin, HISTOGRAM_2D_SPEC* input, int steps, int bunchID);



#endif
