/**
 *
 */


#include "reader2.h"


//
//
//

void init2DHistInput(HISTOGRAM_2D_SPEC* input, int firstCol, int secondCol, const char* n1, const char* n2)
{
  input->_c1 = firstCol, input->_c2 = secondCol; 

  input->_maxHits=0;
  input->_statC1._min = 9999; input->_statC1._max = -9999;
  input->_statC2._min = 9999; input->_statC2._max = -9999;

  sprintf(input->_name1, "%s", n1);
  sprintf(input->_name2, "%s", n2);
}

 
void intervine(int* h1, uint64_t h1size, int* h2, uint64_t h2size, int bin, uint64_t result[bin][bin])
{
  if (h1size != h2size) {
    return;
  }

  uint64_t i,j=0;
  for (i=0; i<h1size; i++) {
      result[h1[i]][h2[i]] ++;
  }
}

void binMap(DATA* v, int* result, int bin)
{
  uint64_t i;
  long_double_t interval = (v->_max - v->_min)/bin;

  for (i=0; i<v->_size; i++) {
    result[i] = (int)((v->_data[i] - v->_min)/interval);    
    
    /*
    if ((i <= 2) || (i >= v->_size -5))  {
      //printf(" .... %dth,  %.8f -> %.8f, vs interval= %.8f, data=%.8f\n", i,  v->_max, v->_min, interval,  v->_data[i]);
      printf("     result[%d] = %d  \n", i, result[i]);
    }
    */
    
  }        

  //printf("\n");
}

void calculate2DHistogram(DATA* v1, DATA* v2, int bin, uint64_t result[bin][bin], int rank)
{

  int  h1[v1->_size], h2[v2->_size];   // values are from [0,bin]

  binMap(v1, h1, bin);
  binMap(v2, h2, bin);
  if (rank == 0) logTimeMillis("       step-binmap\n");    

  intervine(h1, v1->_size, h2, v2->_size,  bin, result);

  //if (rank == 0) logTimeMillis("       intervine  took: \n");    
  
}

//void getAllCol(long_double_t* data, uint64_t nrows, uint64_t ncols, COL_DATA* stat)
/*
void getAllCol(long_double_t* data, uint64_t nrows, uint64_t ncols, COL_DATA* stat)
{
  uint64_t i=0;
  uint64_t r=0;
  for (i=0; i<ncols; i++) {
    //stat[i]._colData=calloc(nrows, sizeof(long_double_t));
    if (stat[i]._colData != NULL) {
      free(stat[i]._colData);
    }
    //stat[i]._colData = getCol(data, nrows, ncols, i, &(stat[i]._colStat));

    stat[i]._colData = calloc(nrows, sizeof(long_double_t));
		   
    stat[i]._colStat._min=9999;
    stat[i]._colStat._max=-9999;
  }
    
  for (r=0; r<nrows; r++) {
      long_double_t v[6] = {data[r*ncols],   data[r*ncols+1], data[r*ncols+2],  
			    data[r*ncols+3], data[r*ncols+4], data[r*ncols+5]};
      for (i=0; i<ncols; i++) {
	stat[i]._colData[r] = v[i];
	stat[i]._colStat._min = MIN(v[i], stat[i]._colStat._min);
	stat[i]._colStat._max = MAX(v[i], stat[i]._colStat._max);
      }
  }
}
*/

void getAllCol(long_double_t* data, uint64_t nrows, uint64_t ncols, COL_DATA* stat)
{
  uint64_t i=0;
  for (i=0; i<ncols; i++) {
    //stat[i]._colData=calloc(nrows, sizeof(long_double_t));
    if (stat[i]._colData != NULL) {
      free(stat[i]._colData);
    }
    stat[i]._colData = getCol(data, nrows, ncols, i, &(stat[i]._colStat));
  }
}



long_double_t* getCol(long_double_t* data, uint64_t nrows, uint64_t ncols, uint64_t whichCol, ARRAY_STAT* stat)
{
  long_double_t* col = calloc(nrows, sizeof(long_double_t));
		   
  stat->_min=9999;
  stat->_max=-9999;

  uint64_t i=0;
  for (i=0; i<nrows; i++) {
    col[i] = data[i*ncols + whichCol];
    stat->_min = MIN(col[i], stat->_min);
    stat->_max = MAX(col[i], stat->_max);
  }

  return col;
}

//void process2D(long_double_t* data, uint64_t nrows, uint64_t ncols, int rank, int nbin, HISTOGRAM_2D_SPEC* input, int steps) {
void process2D(COL_DATA* coldata, uint64_t nrows, uint64_t ncols, int rank, int nbin, HISTOGRAM_2D_SPEC* input, int steps, int bunchID) {

    char filename[100];
    sprintf(filename, "hist2d%d.%s.%s.%d", bunchID, input->_name1, input->_name2, steps);

    if (rank == 0) {
      printf("processing cols: %s, %s with %d bins \n",  input->_name1, input->_name2, nbin);
    }

    //printf(" rank %d, process: %ld, %ld \n", rank, nrows, ncols);
    int i,j;
    
    int c1 = input->_c1;
    int c2 = input->_c2;    
   
#ifdef NEVER   
    printf("    [");
    for (i = 0; i < 6; i++) {
      printf ("%.4f ", data[i]);      
    }
    printf("]\n");
#endif

    FILE* f=fopen(filename, "w");    
    
    //long_double_t *x, *Px, *y, *Py, *z, *Pz=0;

#ifdef NEVER
    long_double_t *coldata1, *coldata2=0;
    ARRAY_STAT statAttr1, statAttr2;
    coldata1 = getCol(data, nrows, ncols, c1, &statAttr1);
    coldata2 = getCol(data, nrows, ncols, c2, &statAttr2);

    long_double_t minAttr1=0, maxAttr1=0; 
    MPI_Allreduce(&(statAttr1._min), &minAttr1, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&(statAttr1._max), &maxAttr1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    long_double_t minAttr2, maxAttr2; 
    MPI_Allreduce(&(statAttr2._min), &minAttr2, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&(statAttr2._max), &maxAttr2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#else
    long_double_t *coldata1 = coldata[c1]._colData;
    long_double_t *coldata2 = coldata[c2]._colData;
    //ARRAY_STAT statAttr1 = coldata[c1]._colStat;
    //ARRAY_STAT statAttr2 = coldata[c2]._colStat;

    long_double_t minAttr1=0, maxAttr1=0; 
    MPI_Allreduce(&(coldata[c1]._colStat._min), &minAttr1, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&(coldata[c1]._colStat._max), &maxAttr1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    long_double_t minAttr2, maxAttr2; 
    MPI_Allreduce(&(coldata[c2]._colStat._min), &minAttr2, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&(coldata[c2]._colStat._max), &maxAttr2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#endif
    input->_statC1._min = MIN(minAttr1,  input->_statC1._min);      input->_statC1._max = MAX(maxAttr1, input->_statC1._max);
    input->_statC2._min = MIN(minAttr2,  input->_statC2._min);      input->_statC2._max = MAX(maxAttr2, input->_statC2._max);


    if (rank == 0) {
      printf("       test:  attr1 min/max: %.8f/%.8f %s\n", minAttr1, maxAttr1, filename);
      printf("       test:  attr2 min/max: %.8f/%.8f %s\n", minAttr2, maxAttr2, filename);
    }

    DATA v1, v2;
    v1._data = coldata1; v1._size = nrows, v1._min = minAttr1, v1._max = maxAttr1;
    v2._data = coldata2; v2._size = nrows, v2._min = minAttr2, v2._max = maxAttr2;

    //if (rank == 0)
    //printf(" rank %d: size of v1/v2=%ld, %ld  bin=%d\n", rank, v1._size, v2._size, nbin);


    uint64_t  hi[nbin][nbin];
    uint64_t  result[nbin][nbin];
    
    for (i=0; i<nbin; i++) {
      for (j=0; j<nbin; j++) {
	result[i][j] = 0;
	hi[i][j]     = 0;
      }
    }

    calculate2DHistogram(&v1, &v2, nbin, result, rank);

    MPI_Reduce(&result[0][0], &hi[0][0], nbin * nbin, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) 
      logTimeMillis("        step MPIReduce");

    if (rank == 0) {
       fprintf(f, "%.8f %.8f ",  minAttr1, maxAttr1);
       fprintf(f, "%.8f %.8f\n", minAttr2, maxAttr2);
       uint64_t counter=0;
       for (i=0; i<nbin; i++) {
	 for (j=0; j<nbin; j++) {
	   counter+= hi[i][j];
	   fprintf(f, "%d %d %ld \n", i, j, hi[i][j]);
	   input->_maxHits = MAX(input->_maxHits, hi[i][j]);
	 }
       }
       printf(" max hits: %ld \n", input->_maxHits);
    }

    //free(coldata1);
    //free(coldata2);

    fclose(f);
}


