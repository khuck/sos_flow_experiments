/*
    originally by Ji Qiang in Fortran
    ported to C by Junmin
    calculate averaged <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
    <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance from
    multiple bunch/bin.
*/

#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "reader2.h"
/*
*/


//#define SQUARE(x) (x*x)
#define SQUARE(x) (pow(x, 2))


//diagnose
void clear(DIAGNOSE* d)
{
  int c,p=0;
  for (c=0; c<NCOL; c++) {
    for (p=0; p<NPOWER; p++) {
      d->_attr[c][p] = 0;
    }
    d->_maxByCol[c] = 0.0;
  }

  d->_xpx = 0.0;
  d->_ypy = 0.0;
  d->_zpz = 0.0;
  d->_gam[0] = 0; d->_gam[1] = 0; // gam & gam2avg
  d->_lrmax=0.0;
  d->_bunchSize=0;    
}



void collect0(long_double_t* all, COL_DATA* coldata, int row)
{
  int col, p;


  for (col=0; col<NCOL; col++) {
    for (p=0; p<NPOWER; p++) {
      if (p == 0) {
	all[col*NPOWER+p] += coldata[col]._colData[row];
      } else if (col < 4)  {
	long_double_t v = coldata[col]._colData[row];
	all[col*NPOWER+p] += pow(v, p+1.0);
      }
      //all[col*NPOWER+p] += pow(coldata[col]._colData[row], p+1);
    }
  }

  all[NCOL*NPOWER  ] += coldata[getPos("x")]._colData[row] * coldata[getPos("px")]._colData[row];
  all[NCOL*NPOWER+1] += coldata[getPos("y")]._colData[row] * coldata[getPos("py")]._colData[row];
#ifdef NEVER
  all[NCOL*NPOWER+2] += coldata[getPos("z")]._colData[row] * coldata[getPos("pz")]._colData[row];
#endif

  all[NCOL*NPOWER+3] += sqrt(1.0 + (SQUARE(coldata[getPos("px")]._colData[row]) +
				    SQUARE(coldata[getPos("py")]._colData[row]) +
				    SQUARE(coldata[getPos("pz")]._colData[row])));

}

void collect(long_double_t* all, COL_DATA* coldata, int row, long_double_t* maxlr, long_double_t maxByCol[6])
{
  long_double_t x  = coldata[0]._colData[row];
  long_double_t px = coldata[1]._colData[row];
  long_double_t y  = coldata[2]._colData[row];
  long_double_t py = coldata[3]._colData[row];
  long_double_t z  = coldata[4]._colData[row];
  long_double_t pz = coldata[5]._colData[row];

#ifdef NEVER
  all[0]+= x;   all[1] += x*x;    all[2] += pow(x,3);   all[3] += pow(x,4);
  all[4]+= px;  all[5] += px*px;  all[6] += pow(px,3);  all[7] += pow(px,4);
  all[8]+= y;   all[9] += y*y;    all[10] += pow(y,3);  all[11] += pow(y,4);
  all[12]+= py; all[13] += py*py; all[14] += pow(py,3); all[15] += pow(py,4);

  all[NCOL*NPOWER  ] += x*px;
  all[NCOL*NPOWER+1] += y*py;

  //all[NCOL*NPOWER+2] += coldata[getPos("z")]._colData[row] * coldata[getPos("pz")]._colData[row];

  all[NCOL*NPOWER+3] += sqrt(1.0 + px*px + py*py + pz*pz);
#else

  long_double_t pow1[4] = {x, px, y, py};
  long_double_t pow2[4] = {x*x, px*px, y*y, py*py};
  long_double_t pow3[4] = {x*pow2[0], px*pow2[1], y*pow2[2], py*pow2[3]};
  long_double_t pow4[4] = {pow2[0]*pow2[0], pow2[1]*pow2[1], pow2[2]*pow2[2], pow2[3]*pow2[3]};


  all[0] += x;   all[1]  += pow2[0];    all[2]  += pow3[0];  all[3]  += pow4[0];
  all[4] += px;  all[5]  += pow2[1];    all[6]  += pow3[1];  all[7]  += pow4[1];
  all[8] += y;   all[9]  += pow2[2];    all[10] += pow3[2];  all[11] += pow4[2];
  all[12]+= py;  all[13] += pow2[3];    all[14] += pow3[3];  all[15] += pow4[3];
  all[16]+= z;
  all[20]+= pz;

  all[NCOL*NPOWER  ] += x*px;
  all[NCOL*NPOWER+1] += y*py;

  //all[NCOL*NPOWER+2] += coldata[getPos("z")]._colData[row] * coldata[getPos("pz")]._colData[row];

  all[NCOL*NPOWER+3] += sqrt(1.0 + px*px + py*py + pz*pz);

  // compute for max
  *maxlr = MAX(pow2[0] + pow2[2], *maxlr);
  int c;
  for (c=0; c<4; c++) { // z pz needs to wait
    maxByCol[c] = MAX(maxByCol[c], fabs(pow1[c]));
  }

#endif
}


long_double_t colQuad(DIAGNOSE* bunchWise, int col)
{

  long_double_t v0  = bunchWise->_attr[col][0];
  long_double_t sq  = bunchWise->_attr[col][1];
  long_double_t cub = bunchWise->_attr[col][2];
  long_double_t fth = bunchWise->_attr[col][3];

  if (col>=4) {
    return sqrt(sqrt(fth));
  }

  long_double_t v04 = sqrt(sqrt(fabs(fth-4*cub*v0 + 6*sq*v0*v0 -3*pow(v0, 4.0))));

  return v04;
}

void getQuads(DIAGNOSE* bunchWise, long_double_t out[6])
{
  int i=0;
  for (i=0; i<6; i++) {
    out[i] = colQuad(bunchWise, i);
  }
}

long_double_t colCube(DIAGNOSE* bunchWise, int col)
{
  long_double_t cub = bunchWise->_attr[col][2];

  if (col >= 4) {
    if (cub > 0) {
      return pow(cub, 1.0/3.0);
    } else {
      return -(pow(fabs(cub), 1.0/3.0));
    }
  }

  long_double_t sq = bunchWise->_attr[col][1];
  long_double_t v0 = bunchWise->_attr[col][0];
  long_double_t v03 = pow(fabs(cub-3*sq*v0+2*pow(v0, 3.0)), (1.0/3.0));

  return v03;
}

void getCubes(DIAGNOSE* bunchWise, long_double_t out[6])
{
  int i=0;
  for (i=0; i<6; i++) {
      out[i] = colCube(bunchWise, i);
      //printf(" cube[%d]= %.10e   ", i, out[i]);
  }
  //printf("\n");
}

void getAlongAxis(DIAGNOSE* bunchWise, int axisPos, long_double_t out[6])
{
  long_double_t sqsum  = bunchWise->_attr[axisPos][1]   - pow(bunchWise->_attr[axisPos][0],2.0);
  long_double_t sqsumP = bunchWise->_attr[axisPos+1][1] - pow(bunchWise->_attr[axisPos+1][0],2.0);  
  long_double_t npn    =-bunchWise->_attr[axisPos][0] * bunchWise->_attr[axisPos+1][0];
  if (axisPos == 0) {
    npn += bunchWise->_xpx;
  } else if (axisPos == 2) {
    npn += bunchWise->_ypy;
  } else {
    npn    = bunchWise->_zpz;
    sqsum  = bunchWise->_attr[axisPos][1];
    sqsumP = bunchWise->_attr[axisPos+1][1];
  }

  long_double_t eps2 = (sqsum*sqsumP-npn*npn);
  long_double_t ep = sqrt(MAX(eps2, 0.0));

  out[0] = bunchWise->_attr[axisPos][0];
  out[1] = sqrt(fabs(sqsum));
  out[2] = bunchWise->_attr[axisPos+1][0];  
  out[3] = sqrt(fabs(sqsumP));
  out[4] = -npn;
  out[5] = ep;

  /*
  int i=0;
  for(i=0; i<6;i++) {
    printf(" out[%d]= %.10e   ", i, out[i]);
  }
  printf("\n");
  */
  //write(24,102)z,z0avg*xl,x0*xl,xrms*xl,px0,pxrms,-xpx*xl,epx*xl    
} 

void printDiagnose(DIAGNOSE* d)
{
  
  int i,c,p=0;
  for (c=0;c<6;c++) {
    for (p=0;p<4; p++) {
      printf("\t%.15E\t",d->_attr[c][p]);
    }
    printf("\n");
  }
  
  printf("\t[xpx]=%.15E, [ypy]=%.15E, [zpz]= %.15E gamma={%.15E, %.15E} \n", d->_xpx, d->_ypy, d->_zpz, d->_gam[0], d->_gam[1]);	 
  for (i=0; i<6; i++) {
    printf("\t... [%.15E]\t", d->_maxByCol[i]);
  }
  printf("\n");
  printf("\t.... bunch size=%d\n\n", d->_bunchSize);
  
}

void finishDiagnose(DIAGNOSE* diagOfAllBunchs, int numBunch)
{
  int i,c,p=0;
  /*
  for (i=0; i<numBunch; i++) {
    printf("BUNCH: %d === \n", i);
    printDiagnose(&(diagOfAllBunchs[i]));
  }
  */
  
  long_double_t   in_dt=4.0e-12; // input File
  long_double_t   in_mass = 0.511005e06; // input File
  uint64_t in_speedLight=299792458; // given
  long_double_t in_rad2deg = 180/(2*asin(1.0));

  long_double_t qmc =  in_mass/1.0e6;
  long_double_t xl = (double)in_speedLight*in_dt;
  long_double_t xt = in_rad2deg;

  //printf("... xl = %.10e \n", xl);

  long_double_t z = 0.0; // related to timestep
  
  // now calculate
  long_double_t gamavg=0;
  uint64_t nptot = 0;
  long_double_t glrmax = 0;
  long_double_t gam2avg = 0;
  //long_double_t x0 = 0;

  DIAGNOSE bunchWise;
  clear(&bunchWise);

  for (i=0; i<numBunch; i++) {
    for (c=0;c<6;c++) {
      for (p=0;p<4; p++) {	
	bunchWise._attr[c][p] += diagOfAllBunchs[i]._attr[c][p];
      }
      bunchWise._maxByCol[c] = MAX(diagOfAllBunchs[i]._maxByCol[c], bunchWise._maxByCol[c]);
    }

    bunchWise._xpx += diagOfAllBunchs[i]._xpx;
    bunchWise._ypy += diagOfAllBunchs[i]._ypy;
    bunchWise._zpz += diagOfAllBunchs[i]._zpz;
    bunchWise._gam[0] += diagOfAllBunchs[i]._gam[0];
    bunchWise._gam[1] += diagOfAllBunchs[i]._gam[1];

    bunchWise._bunchSize += diagOfAllBunchs[i]._bunchSize;
    bunchWise._lrmax = MAX(glrmax, diagOfAllBunchs[i]._lrmax);
  }
 
  
  //printf("====== check sum of BUNCHs=== \n");
  //printDiagnose(&bunchWise);


  //withAllOthers(&bunchWise);

  for (c=0;c<6;c++) {
    for (p=0;p<4; p++) {
      bunchWise._attr[c][p] = bunchWise._attr[c][p]/(double)bunchWise._bunchSize;
    } 
  }   
  bunchWise._gam[0] = bunchWise._gam[0]/(double)bunchWise._bunchSize;
  bunchWise._gam[1] = bunchWise._gam[1]/(double)bunchWise._bunchSize;
  
  bunchWise._xpx = bunchWise._xpx/(double)bunchWise._bunchSize;
  bunchWise._ypy = bunchWise._ypy/(double)bunchWise._bunchSize;
  bunchWise._zpz = bunchWise._zpz/(double)bunchWise._bunchSize;
  gam2avg = gam2avg/(double)bunchWise._bunchSize;


  //printf("BUNCHWISE=== \n");
  //printDiagnose(&bunchWise);

  gamavg = bunchWise._gam[0];
  long_double_t energy = qmc*(gamavg-1.0f);
  long_double_t bet = sqrt(1.0f-(1.0f/(gamavg*gamavg)));  

  gam2avg = bunchWise._gam[1];
  long_double_t gamdel = sqrt(gam2avg);
 
  printf("fort 18: %.15E, %.15E, %.15E, %.15E, %.15E, %.15E, %.15E\n", z, bunchWise._attr[4][0]*xl, gamavg, energy, bet, sqrt(bunchWise._lrmax)*xl, gamdel);
  //write(18,100)z,z0avg*xl,gam,energy,bet,sqrt(glrmax)*xl,gamdel

  long_double_t axisOut[6];
  getAlongAxis(&bunchWise, 0, axisOut);
  printf("fort 24: %.15E, %.15E, %.15E, %.15E, %.15E, %.15E, %.15E, %.15E\n", 
	 z, bunchWise._attr[4][0]*xl, axisOut[0]*xl, axisOut[1]*xl, axisOut[2], axisOut[3], axisOut[4]*xl, axisOut[5]*xl);

  getAlongAxis(&bunchWise, 2, axisOut);
  printf("fort 25: %.15E, %.15E, %.15E, %.15E, %.15E, %.15E, %.15E, %.15E\n", 
	 z, bunchWise._attr[4][0]*xl, axisOut[0]*xl, axisOut[1]*xl, axisOut[2], axisOut[3], axisOut[4]*xl, axisOut[5]*xl);

  getAlongAxis(&bunchWise, 4, axisOut);
  printf("fort 26: %.15E, %.15E, %.15E, %.15E, %.15E, %.15E, %.15E\n", 
	 z, bunchWise._attr[4][0]*xl, /*axisOut[0]*xl,*/ axisOut[1]*xl, axisOut[2], axisOut[3], axisOut[4]*xl, axisOut[5]*xl);

  //  write(24,102)z,z0avg*xl,x0*xl,xrms*xl,px0,pxrms,-xpx*xl,epx*xl
  //  write(25,102)z,z0avg*xl,y0*xl,yrms*xl,py0,pyrms,-ypy*xl,epy*xl
  //  write(26,100)z,z0avg*xl,      zrms*xl,pz0,pzrms,-zpz*xl,epz*xl

  printf("fort 27: %.15E, %.15E, %.15E, %.15E, %.15E, %.15E, %.15E %.15E\n", z, bunchWise._attr[4][0]*xl, 
	 bunchWise._maxByCol[0]*xl, bunchWise._maxByCol[1],
	 bunchWise._maxByCol[2]*xl, bunchWise._maxByCol[3], bunchWise._maxByCol[4]*xl, bunchWise._maxByCol[5]);
  //write(27,102)z,z0avg*xl,glmax(1)*xl,glmax(2),glmax(3)*xl,&glmax(4),glmax(5)*xl,glmax(6)
    
  // has to be done at server!!
  //write(28,101)z,z0avg*xl,npctmin,npctmax,nptot  

  long_double_t cubes[6];
  getCubes(&bunchWise, cubes);

  printf("fort 29: %.15E, %.15E, %.15E, %.15E, %.15E, %.15E, %.15E, %.15E\n", 
	 z, bunchWise._attr[4][0]*xl, cubes[0]*xl, cubes[1], cubes[2]*xl, cubes[3], cubes[4]*xl, cubes[5]);
  //write(29,102)z,z0avg*xl,x03*xl,px03,y03*xl,py03,z03*xl,&pz03
  
  long_double_t quads[6];
  getQuads(&bunchWise, quads); 

  printf("fort 30: %.15E, %.15E, %.15E, %.15E, %.15E, %.15E, %.15E, %.15E\n\n", 
	 z, bunchWise._attr[4][0]*xl, quads[0]*xl, quads[1], quads[2]*xl, quads[3], quads[4]*xl, quads[5]);
  //write(30,102)z,z0avg*xl,x04*xl,px04,y04*xl,py04,z04*xl,&pz04
}
 
void update(DIAGNOSE* d, long_double_t* raw)
{
  int c,p=0;
  for (c=0; c<NCOL; c++) {
    for (p=0; p<NPOWER; p++) {
      d->_attr[c][p] = raw[NPOWER*c+p];
    }
  }

  d->_xpx = raw[NPOWER*NCOL];
  d->_ypy = raw[NPOWER*NCOL+1];
  d->_zpz = raw[NPOWER*NCOL+2];
  d->_gam[0] = raw[NPOWER*NCOL+3];

  d->_lrmax = raw[NPOWER*NCOL+4];

}

long_double_t getSum(COL_DATA* coldata, uint64_t size, int col)
{
  int i=0;
  long_double_t sum = 0;
  for (i=0; i<size; i++) {
    sum += coldata[col]._colData[i];
  } 
}

void prepareDiagnose(COL_DATA* coldata, uint64_t size, int ncols,  int rank, DIAGNOSE* stat)
{
  int i=0, c=0;
  // save sum of power(n), n=1,4 for each col
  // save sum of xpx/ypy/zpz 
  long_double_t all[NCOL*NPOWER+5] = {0};


  long_double_t maxlr = 0;
  long_double_t maxByCol[NCOL] = {-9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -9999.0};

 
  for (i=0; i<size; i++) {
    //collect(all, coldata, i);
    collect(all, coldata, i, &maxlr, maxByCol);    
  }
  /*
  printf(" this run: maxlr=%.10e, maxByCol={%10.e, %10.e, %10.e, %10.e, }\n", maxlr, maxByCol[0], maxByCol[1], maxByCol[3], maxByCol[3]);

     maxlr=0;
     for(c=0; c<4; c++) {
       maxByCol[c]=-9999.0;
     }
  
  for (i=0; i<size; i++) {
    long_double_t lr = SQUARE(coldata[getPos("x")]._colData[i]) + SQUARE(coldata[getPos("y")]._colData[i]);
    maxlr = MAX(maxlr, lr);

    for (c=0; c<4; c++) { // z pz needs to wait
      maxByCol[c] = MAX(maxByCol[c], fabs(coldata[c]._colData[i]));
    }
  }
  
  printf(" double check run: maxlr=%.10e, maxByCol={%10.e, %10.e, %10.e, %10.e, }\n", maxlr, maxByCol[0], maxByCol[1], maxByCol[3], maxByCol[3]);
  */
  /*
  printf("... rank=%d     all[0]  =%.10e\n ", rank, all[0]);
  long_double_t hi=0,hi2=0;
  MPI_Allreduce(&all[0], &hi,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&all[0], &hi2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  printf("... rank=%d, hi=%.10e  hi2=%.10e\n", rank, hi, hi2);
  */

  long_double_t sumOver[NCOL*NPOWER+5]= {0};
  
  //printf("... rank=%d  now all[0]  =%.10e \n", rank, all[0]);
  MPI_Allreduce(&all[0], &sumOver[0],  NCOL*NPOWER+4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //printf("... rank=%d  now sum[0]  =%.10e", rank, sumOver[0]);
  long_double_t glrmax;
  MPI_Allreduce(&maxlr, &glrmax, 1, MPI_DOUBLE, MPI_MAX,  MPI_COMM_WORLD);

  sumOver[NCOL*NPOWER+4] = glrmax;

  update(stat, sumOver);

  MPI_Allreduce(&maxByCol[0], &(stat->_maxByCol[0]), 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  //printDiagnose(stat);
}

