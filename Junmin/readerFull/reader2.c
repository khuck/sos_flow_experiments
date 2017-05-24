/**
 *
 */

#include "reader2.h"

const char ATTR_NAMES[6] [20] = {"x", "Px", "y", "Py", "z", "Pz"};
uint64_t _lastMeasuredMillis = (uint64_t)0;

uint64_t getCurrentTimeMillis() {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return (tv.tv_sec*(uint64_t)1000000+tv.tv_usec)/(uint64_t)1000;
}


//
//
//
int getPos(const char* colName) {
  if ((strcmp(colName, "x") == 0) || (strcmp(colName, "X") == 0)) {
    return 0;
  }
 
  if ((strcmp(colName, "Px") == 0) || (strcmp(colName, "PX") == 0) || (strcmp(colName, "px") == 0) ) {  
    return 1;
  }

  if ((strcmp(colName, "y") == 0) || (strcmp(colName, "Y") == 0)) {
    return 2;
  }

  if ((strcmp(colName, "Py") == 0) || (strcmp(colName, "PY") == 0) || (strcmp(colName, "py") == 0) ) {  
    return 3;
  } 

  if ((strcmp(colName, "z") == 0) || (strcmp(colName, "Z") == 0)) {
    return 4;
  }

  if ((strcmp(colName, "Pz") == 0) || (strcmp(colName, "PZ") == 0) || (strcmp(colName, "pz") == 0) ) {  
    return 5;
  }

  printf("Error: Invalid colName %s \n", colName);
  return -1;
}


void printMe(int rank, HISTOGRAM_2D_SPEC* input /*, char* name*/)
{
  if (rank == 0) {
    printf("Summary %s/%s: %.3f/%.3f, %.3f/%.3f maxHits=%ld\n", input->_name1, input->_name2, input->_statC1._min, input->_statC1._max, input->_statC2._min, input->_statC2._max, input->_maxHits);
  }

}

void findAttrValue(ADIOS_FILE* adios_handle, const char* attrName, int rank, int* attrValue)
{
  ADIOS_VARINFO* avi = adios_inq_var (adios_handle, attrName);
  if (!avi){
    printf("rank %d: Quitting ... (%d) %s\n", rank, adios_errno, adios_errmsg());
    adios_read_close(adios_handle);                                                               
    return;
  }
  *attrValue = *((int*)avi->value);
  
  /*
  if (rank == 0)
    printf("%s=%ld\n", attrName, *attrValue);
  */

  adios_free_varinfo(avi);
  avi = NULL;
}

uint64_t getData(long_double_t** t, ADIOS_FILE* adios_handle, int numCores, int rank, int NX, int bunch, uint64_t* bunchSize)
{
    char varName[40];

    if (bunch < 10) {
      sprintf(varName, "particles0000%d", bunch);
    } else if (bunch < 100) {
      sprintf(varName, "particles000%d",  bunch);
    } else if (bunch < 1000) {
      sprintf(varName, "particles00%d",   bunch);
    } else if (bunch < 10000) {
      sprintf(varName, "particles0%d",    bunch);
    } else if (bunch < 100000) {
      sprintf(varName, "particles%d",   bunch);
    } else {
      printf("## Error: Unable to handle large bunch! limit is 100000 ##\n");
      return 0;
    }

    ADIOS_SELECTION *sel = NULL;
    ADIOS_VARINFO* v = adios_inq_var (adios_handle, varName);
    //adios_inq_var_blockinfo (adios_handle, v);
    
    *bunchSize = v->dims[0];

    int sliceBase = numCores;
    uint64_t slice_size = v->dims[0]/(uint64_t)sliceBase;
    
    uint64_t start[2] = { rank * slice_size, 0};
    
    if (rank == sliceBase-1)
      slice_size = slice_size + v->dims[0]%sliceBase;
    
    uint64_t count[2] = {slice_size, NX};
    
    
    printf("....rank %d, start: %d, count=%d , slice=%d\n", rank, start[0], count[0], slice_size);
    sel = adios_selection_boundingbox(2,start, count);
    if( !sel ){
      printf("rank %d: Quitting ... (%d) %s\n", rank, adios_errno, adios_errmsg());
      adios_read_close(adios_handle);                                                               
      return 0;
    }
    //*t = calloc(slice_size*NX, sizeof(long_double_t));
    *t = calloc(slice_size*NX, adios_type_size(v->type, v->value));


    //MPI_Barrier (MPI_COMM_WORLD);

    if (rank == 0) {
        printf ("--------- Step: %d --------------------------------type size=%d\n", adios_handle->current_step, adios_type_size(v->type, v->value));
	//printf("dimension of particle00001: [%d, %d]\n", v->dims[0], v->dims[1]);
	printf("dimension of %s: [%d, %d]\n", varName, v->dims[0], v->dims[1]);
        printf ("..... last: %d, is streaming? %d\n", adios_handle->last_step, adios_handle->is_streaming);
    }

    uint64_t t0 = getCurrentTimeMillis();
    //printf(" RANK=%d, step=%d ++\n", rank, steps);
    if (adios_schedule_read (adios_handle, sel, varName, 0, 1, *t) != 0) {
      printf("rank %d: Quitting ...(%d) %s\n", rank, adios_errno, adios_errmsg());
      adios_read_close(adios_handle);
      return 0;
    }

    adios_perform_reads (adios_handle, 1);


    uint64_t t1 = getCurrentTimeMillis();
    uint64_t diff=t1-t0;
    uint64_t maxRead =0 ;
    MPI_Reduce(&diff, &maxRead, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
      printf("FYI reading took %llu millis, max=%llu\n", t1-t0, maxRead);
    }
    
    adios_selection_delete(sel); sel=NULL;
    adios_free_varinfo(v);
    uint64_t t2 = getCurrentTimeMillis();
    adios_release_step(adios_handle); // added via norbert's suggestion
    uint64_t t3 = getCurrentTimeMillis();
    if (rank == 0) {
      printf("FYI reading took %llu millis, max=%llu, t2-t1=%llu millis, t3-t2=%llu secs\n", t1-t0, maxRead, t3-t2, t2-t1);
    }
  return slice_size;
}

void printAllCol(COL_DATA* allColData, uint64_t nrows, int ncols)
{
  int c, r;
  for (r=0; r<nrows; r++) {
    for (c=0; c<ncols; c++) {
      printf("%.10e ", allColData[c]._colData[r]);
      //printf(" %Lf ", allColData[c]._colData[r]);
      //printf("*%.10llf ", allColData[c]._colData[r]);
    }
    printf("\n");
  }   
}

void freeCol(COL_DATA* allColData, int whichCol)
{
  if (allColData[whichCol]._colData != NULL) {
    free(allColData[whichCol]._colData);
  } 
}

void freeAllCol(COL_DATA* allColData, int ncols)
{
  int c;
  for (c=0; c<ncols; c++) {
    freeCol(allColData, c);
  }
}

long_double_t getDistance(long_double_t a, long_double_t b, long_double_t c)
{
  return (sqrt(1+a*a + b*b + c*c));
} 


long_double_t kahanSum(long_double_t input[], uint64_t size)
{
  uint64_t i=0;
  long_double_t c = 0.0;
  long_double_t sum = input[0];
  for (i=1; i<size; i++) {
    long_double_t y = input[i] - c;
    long_double_t t = sum + y;
    c = (t-sum) - y;
    sum = t;
  }
  printf(".......kahanSum=%.20E\n", sum);
  return sum;
}

void sumPowers(long_double_t p1, long_double_t store[])
{
  long_double_t p2 = p1*p1;
  long_double_t p3 = fabs(p1*p2);
  long_double_t p4 = p2*p2;
  store[0]+=p2; store[1]+=p3; store[2] +=p4;

}


int readFile(char* fileName, int rank, int numCores, ADIOS_FILE* adios_handle, int nbin,   HISTOGRAM_2D_SPEC inputColumns[], int num2DPairs, int numBunch)
{
  //int gysize, NX, NY, GX, GY, OX, OY;
  int gysize, NX,  GY;

  long_double_t* rawData = NULL;
  if ( !adios_handle){
    printf("Unable to read file \n");
    return 0;
  }

  if (rank == 0) 
    printf ("%s... file is read. numCores=%d\n", fileName, numCores);

  // define portions of data how they will be read
  ADIOS_VARINFO *avi = NULL;

  // for storing the variables
  
  int step = 0;
  
  findAttrValue(adios_handle, "GY00001", rank, &gysize);

  // if I run the more readers than writers; just release
  // the excessive readers

  findAttrValue(adios_handle, "NX00001", rank, &NX);
  //findAttrValue(adios_handle, "OX00001", rank, &OX);
  //findAttrValue(adios_handle, "OY00001", rank, &OY);


  if (rank == 0) {
    printf(" rank=%d, NX=%d, GY=%d\n", rank, NX, gysize);
  }
  

  logTimeMillis("      starting to read.\n");
  int steps=0;
  int i;
  float timeout_sec=1.0;


  DIAGNOSE stat[numBunch];
  COL_DATA zCol[numBunch];
  COL_DATA pzCol[numBunch];
  COL_DATA gmCol[numBunch];
  uint64_t sizeAtThisCore[numBunch];
  
  uint64_t t00=0, t0=0, t1=0, t2=0;
  if (rank == 0) logTimeMillis(" step 0\n");
  while (adios_errno != err_end_of_stream) {
    uint64_t stageStart = _lastMeasuredMillis;
    //if (_lastMeasuredMillis == 0) {
      
    steps++; // steps start counting from 1
    int i=0;

    int bunchID=1;    
    for (bunchID=0; bunchID<numBunch; bunchID++) {
         uint64_t bunchSize = 0; // at this step
	 if (rawData != NULL) {
	   free(rawData); rawData = NULL;
	 }
	 uint64_t slice_size  = getData(&rawData, adios_handle, numCores, rank, NX, bunchID+1, &bunchSize);	 
	 if (rank == 0) logTimeMillis(" step 1-read.raw\n");
	 sizeAtThisCore[bunchID] = slice_size;

	 clear(&stat[bunchID]);
	 stat[bunchID]._bunchID   = bunchID+1;
	 stat[bunchID]._bunchSize = bunchSize;
	 stat[bunchID]._step     = steps;

	 t00 = getCurrentTimeMillis();
	 COL_DATA allColData[NX];
	 for (i=0; i<NX; i++) {
	   allColData[i]._colData = NULL;
	 }	 
	 getAllCol(rawData, slice_size, NX, allColData);
	 if (rank == 0) {
	   //printf("  t00 = %ld vs %ld\n", t00, _lastMeasuredMillis);
	   logTimeMillis(" step 2-getCol\n");
	 }
	 //printAllCol(allColData, slice_size, NX);
	 
	  MPI_Barrier (MPI_COMM_WORLD);

#ifdef NEVER // no hist computation now
	  i=0;
	  while (i < num2DPairs) {
	    process2D(allColData, slice_size, NX, rank, nbin, &(inputColumns[i]), steps, bunchID+1);
	   i++;
	 }	 
	  if (rank == 0) logTimeMillis(" step 2-2d\n");
#endif
#ifndef DIAG
	  t0 = getCurrentTimeMillis();
	  //if (rank == 0) printf("   t0 = %ld, vs %ld \n", t0, _lastMeasuredMillis);
	  // prepare all the diagnose stat for this step	 
	  prepareDiagnose(allColData, slice_size, NX, rank, &stat[bunchID]);
	  t1 = getCurrentTimeMillis();	  
	  if (rank == 0) logTimeMillis(" step 3-pre\n");
#endif
	 //freeAllCol(allColData, NX);
	 freeCol(allColData, 0);freeCol(allColData, 2);
	 
	 gmCol[bunchID]._colData = calloc(slice_size, sizeof(long_double_t));
	 for (i=0; i<slice_size; i++) {
	   gmCol[bunchID]._colData[i] = getDistance(allColData[1]._colData[i], allColData[3]._colData[i], allColData[5]._colData[i]);
	 }

	 freeCol(allColData, 1);freeCol(allColData, 3);
	 zCol[bunchID]._colData  = allColData[4]._colData;
	 pzCol[bunchID]._colData = allColData[5]._colData;
	 t2 = getCurrentTimeMillis();	  
	 if (rank == 0) logTimeMillis(" step 4-diag\n");
    }

#ifndef DIAG
    //printDiagnose(&(stat[0]));
    // adjustment for getting effective scientific digits
    long_double_t z0avg=0.0, pz0avg=0.0, gam0avg=0.0;
    uint64_t sumPtl=0;
    for (bunchID=0; bunchID<numBunch; bunchID++) {
       z0avg  += stat[bunchID]._attr[4][0] ;
       pz0avg += stat[bunchID]._attr[5][0];
       gam0avg += stat[bunchID]._gam[0];
       sumPtl += stat[bunchID]._bunchSize;
    }

    printf(" checking z0avg: %.20E/%llu \n", z0avg, sumPtl);
    z0avg = z0avg/(double)sumPtl; pz0avg = pz0avg/(double)sumPtl; gam0avg=gam0avg/(double)sumPtl;
    
    long_double_t maxZ=-9999.0, maxPz=-9999.0;
    for (bunchID=0; bunchID<numBunch; bunchID++) {
         long_double_t zsq[3]={0}, pzsq[3]={0};
	 long_double_t zpz  = 0, gam2=0;
	 int k=0;
	 for (i=0; i<sizeAtThisCore[bunchID]; i++) {
	      long_double_t z1 = zCol[bunchID]._colData[i] - z0avg;
	      sumPowers(z1, zsq);
	      long_double_t pz1 = pzCol[bunchID]._colData[i] - pz0avg;
	      sumPowers(pz1, pzsq);

	      maxZ   = MAX(maxZ, fabs(z1));
	      maxPz  = MAX(maxPz, fabs(pz1));
	      zpz   += z1*pz1;

	      //maxZ  = MAX(maxZ,  fabs(zCol[bunchID]._colData[i]  - z0avg));
	      //maxPz = MAX(maxPz, fabs(pzCol[bunchID]._colData[i] - pz0avg));
	      //zpz  += (zCol[bunchID]._colData[i] - z0avg) * (pzCol[bunchID]._colData[i] - pz0avg);
	      gam2 += (gmCol[bunchID]._colData[i] - gam0avg) * (gmCol[bunchID]._colData[i] - gam0avg); 
	 }
	 
	 //printf("rank=%d p2=%.20E p3=%20E p4=%.20E\n", rank, zsq[0], zsq[1], zsq[2]);
	 
	 long_double_t sumzsq[3] = {0};
	 long_double_t sumpzsq[3] = {0};
	 long_double_t sumgam2, sumzpz, gMaxZ, gMaxPz;

	 MPI_Reduce(&zsq,  &sumzsq,  3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	 MPI_Reduce(&pzsq, &sumpzsq, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	 MPI_Reduce(&gam2, &sumgam2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	 MPI_Reduce(&zpz,  &sumzpz,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 

	 MPI_Reduce(&maxZ,  &gMaxZ,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 
	 MPI_Reduce(&maxPz, &gMaxPz, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 

	 for (k=2; k<5; k++) {
	   stat[bunchID]._attr[4][k-1] = sumzsq[k-2];
	   stat[bunchID]._attr[5][k-1] = sumpzsq[k-2];
	 }       

	 stat[bunchID]._zpz = sumzpz;
	 stat[bunchID]._gam[1] = sumgam2;
	 stat[bunchID]._maxByCol[4] = gMaxZ;
	 stat[bunchID]._maxByCol[5] = gMaxPz;
    }

    if (rank == 0) logTimeMillis(" step 5\n");
    uint64_t t3 = getCurrentTimeMillis();	  

    if (rank == 0) {
      finishDiagnose(stat, numBunch);
    }
    //logTimeMillis(" step 6\n");
    uint64_t t4 = getCurrentTimeMillis();	  
    uint64_t maxDiag=0;
    uint64_t diff = t4-t0;
    MPI_Reduce(&diff, &maxDiag, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      printf("  t0=%llu, t1 = %llu mill, t2 = %llu mill, t3 = %llu mill, t4 = %llu mill, total=%llu mill, max=%llu\n", t0-t00, t1-t0, t2-t1, t3-t2, t4-t3, t4-t0, maxDiag);
    }
    //logTimeMillis(" step 7\n");
#endif    

    //printf("\n free auxillary columns!!! \n");
    for (i=0; i<numBunch; i++) {
      free(zCol[i]._colData);free(pzCol[i]._colData);free(gmCol[i]._colData);
    }
    if (rank == 0) logTimeMillis(" step 8\n");
    adios_advance_step (adios_handle, 0, timeout_sec);
    if (rank == 0) logTimeMillis(" step 9\n");
    if (adios_errno == err_step_notready) {      
        printf ("rank %d: No new step arrived within the timeout. Quit. %s\n", rank, adios_errmsg());		
	break; // quit while loop
    }    

    
    if (rank == 0) {
        uint64_t stageEnd = getCurrentTimeMillis();
        printf(" .... so ts=%d, TOOK : %d milli sec \n",  steps, (stageEnd-stageStart)); 
    }  
  }

  printMe(rank, inputColumns);

  //adios_free_varinfo(v);
  //adios_selection_delete(sel); sel=NULL;
  free(rawData); rawData = NULL;
  MPI_Barrier (MPI_COMM_WORLD);
  adios_read_close(adios_handle);


  return 0;	
}


int main (int argc, char **argv) 
{
  int rank =0, size =0;
    MPI_Comm comm = MPI_COMM_WORLD;

    // adios read initialization
    MPI_Init( &argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    int nbin=5;
    if (argc < 4) {
      printf("Usage: %s bin colname1a colname1b ...\n", argv[0]);
      MPI_Finalize();
      return 0;
    } 

    if (argc % 2 > 0) {
      printf("Error: Column names should be in pairs.\n");
      printf("Usage: %s bin colname1a colname1b ...\n", argv[0]);
      MPI_Finalize();
      return 0;
    } 
      

    nbin = atoi(argv[1]);
    
    enum ADIOS_READ_METHOD flexMethod = ADIOS_READ_METHOD_FLEXPATH; 
    int err = adios_read_init_method(flexMethod, comm, "");
    if (err != 0) {
      printf("Unable to set flexpath method \n");
      MPI_Finalize();
      return 0;
    }
    
    _lastMeasuredMillis = getCurrentTimeMillis();
    
    char* startFileName="start00001.bp";
    ADIOS_FILE *adios_handle = adios_read_open(startFileName, flexMethod, comm, ADIOS_LOCKMODE_NONE, 0.0);
    
    logTimeMillis("         file opened. \n ");
    
    int numBunch;
    findAttrValue(adios_handle, "NBunch", rank, &numBunch);
    if (rank == 0) {
      printf(" num of bunch=%d\n", numBunch);
    }

    int total2DPairs= (argc-2)/2;
    HISTOGRAM_2D_SPEC inputColumns[total2DPairs];

    int pos=0;
    while ((pos+1)*2 < argc) {
      char* firstCol  = argv[(pos+1)*2] ;
      char* secondCol = argv[(pos+1)*2+1];

      int c1 = getPos(firstCol);
      int c2 = getPos(secondCol);

      init2DHistInput(&inputColumns[pos], c1, c2, ATTR_NAMES[c1], ATTR_NAMES[c2]);      
      if ((c1 < 0) || (c2 < 0)) {
	MPI_Finalize();
	return 0;
      }
      pos +=1;
    }
    
    readFile(startFileName, rank, size, adios_handle, nbin, inputColumns, total2DPairs, numBunch);
    logTimeMillis("      finished.\n");
    adios_read_finalize_method(flexMethod);                                             

    MPI_Finalize();
}

