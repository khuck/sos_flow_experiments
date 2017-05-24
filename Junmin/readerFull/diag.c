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



int DefColumnX=0;
int DefColumnPx=1;
int DefColumnY=2;
int DefColumnPy=3;
int DefColumnZ=4;
int DefColumnPz=5;

int DefNumCols=6;

struct _colStat 
{
  // note: sum[i] = pow(i+1)
  double  _sum[4];   
};

typedef struct _colStat ColStat;

struct _pointStat 
{
  double _xPx;
  double _yPy;
  double _zPz;
  
  double _max[6]; // localmax
  double _lcrmax;

  double _gamlc;
  double _gam2lc;
};

typedef struct _pointStat PointStat;

//


#define CUBE(x) (x*x*x)
#define SQUARE(x) (x*x)
#define SQUARE_SUM(x, y, z) (x*x + y*y + z*z)
#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)

void init(ColStat* s)
{
  int i=0;
  for (i=0; i<4; i++) {
    s->_sum[i] = 0.0;
  }
}

void SummarizeCol(ColStat* stat, ColStat* sumStat)
{
  MPI_Reduce(stat->_sum, sumStat->_sum, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

double get03(ColStat* sumStat, double den1)
{
  double v0   = (sumStat->_sum[0]) * den1;
  double sqv  = (sumStat->_sum[1]) * den1;
  double cubv = (sumStat->_sum[2]) * den1;

  double base = abs(cubv - 3*sqv*v0 + 2*CUBE(v0));
  double v03 = pow(base, 1.0/3.0);
  return v03;
}

double get04(ColStat* sumStat, double den1)
{
  double v0   = (sumStat->_sum[0]) * den1;
  double sqv  = (sumStat->_sum[1]) * den1;
  double cubv = (sumStat->_sum[2]) * den1;
  double fthv = (sumStat->_sum[3]) * den1;

  return sqrt(sqrt(abs(fthv-4*cubv*v0+6*sqv*v0*v0-3*v0*v0*v0*v0)));
}


void initPoint(PointStat* s)
{
  s->_xPx = 0.0;
  s->_yPy = 0.0;
  s->_zPz = 0.0;

  s->_lcrmax = 0.0;

  int i=0;
  for (i=0; i<6; i++) {
    s->_max[i] = 0.0;
  }
}


double pSquare(double* data, uint64_t pos)
{	    
  return SQUARE(data[pos+DefColumnPx]) + SQUARE(data[pos+DefColumnPy]) + SQUARE(data[pos+DefColumnPz]);
}

void sumtmp3lc(double* data, uint64_t nrows,  double result[3])
{
  uint64_t i=0; 
  
  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;
    
  for (i=0; i<nrows; i++) {
    result[0] += data[i* DefNumCols + DefColumnZ];
    result[1] += data[i* DefNumCols + DefColumnPz];
    result[2] += sqrt(1.0 + pSquare(data, i * DefNumCols));
  }
}

void onPoint(double* data, uint64_t whichrow, PointStat* pointStat, double z0avg, double pz0avg, double gamavg)
{
  double x = data[whichrow * DefNumCols + DefColumnX];
  double px= data[whichrow * DefNumCols + DefColumnPx];

  double y = data[whichrow * DefNumCols + DefColumnY];
  double py= data[whichrow * DefNumCols + DefColumnPy];

  double z = data[whichrow * DefNumCols + DefColumnZ];
  double pz= data[whichrow * DefNumCols + DefColumnPz];

  pointStat->_xPx +=  x * px;
  pointStat->_yPy +=  y * py;
  pointStat->_zPz +=  (z-z0avg)*(pz-pz0avg);

  int i=0;

  double center[6] = {0.0, 0.0, 0.0, 0.0, z0avg, pz0avg};

  for (i=0; i<6; i++) {
    pointStat->_max[i] = MAX(pointStat->_max[i], abs(data[whichrow * DefNumCols + i] - center[i]));
  }

  pointStat->_lcrmax = MAX(pointStat->_lcrmax, x*x + y*y);

  double tmpgam = sqrt(1.0 + SQUARE_SUM(px, py, pz));
  pointStat->_gamlc  += tmpgam;
  
  pointStat->_gam2lc += SQUARE(tmpgam-gamavg);

}

void onColumn(double* data, uint64_t whichrow,  int whichCol, ColStat* stat, double center)
{
    double v = data[whichrow * DefNumCols + whichCol];

    int i=0;
    for (i=0; i<4; i++) {
      if (i > 1) {
	stat->_sum[i] += pow(v - center, i);
      } else {
	stat->_sum[i] += pow(v, i);
      }
    }
}

void SummarizePoint(PointStat* stat, PointStat* sumStat)
{
  double input[5] = {stat->_xPx, stat->_yPy, stat->_zPz, stat->_gamlc, stat->_gam2lc};
  double output[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  MPI_Reduce(input, output, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  sumStat->_xPx = output[0];
  sumStat->_yPy = output[1];
  sumStat->_zPz = output[2];
  sumStat->_gamlc  = output[3];
  sumStat->_gam2lc = output[4];
}


void write18(PointStat* sumPointStat, double den1, double qmc)
{
  double gam = (sumPointStat->_gamlc) * den1;
  double energy = qmc * (gam-1.0);

  double bet = sqrt(1.0 - SQUARE(1.0/gam));

  double gam2avg = (sumPointStat->_gam2lc) * den1;
  double gamdel = sqrt(gam2avg);

  // write (z, z0avg*DefConstXl, gam, energy, bet, sqrt((pointStat._lcrmax))*DefConstXl, gamdel)
}

void write24(ColStat* sumStatX, ColStat* sumStatPx, PointStat* sumPointStat, double den1)
{
  double x0 = (sumStatX->_sum[0]) * den1;
  double sqx = (sumStatX->_sum[1]) * den1;
  double sqsum1 = sqx - x0*x0;

  double px0 = (sumStatPx->_sum[0]) * den1;
  double sqpx = (sumStatPx->_sum[1])* den1;
  double sqsum2 = sqpx - px0*px0;
  
  double xrms = sqrt(abs(sqsum1));
  double pxrms = sqrt(abs(sqsum2));

  double xpx = (sumPointStat->_xPx)* den1 - x0 * px0;

  // write (z, z0avg * DefConstXl, x0*DefConstXl, xrms*DefConstXl, px0, pxrms, -xpx*DefConstXl, epx*DefConstXl);
}

void write25(ColStat* sumStatY, ColStat* sumStatPy, PointStat* sumPointStat, double den1)
{
  double y0  = (sumStatY->_sum[0]) * den1;
  double sqy = (sumStatY->_sum[1]) * den1;
  double sqsum3 = sqy - y0*y0;

  double py0 = (sumStatPy->_sum[0]) * den1;
  double sqpy = (sumStatPy->_sum[1])* den1;
  double sqsum4 = sqpy - py0*py0;
  
  double yrms = sqrt(abs(sqsum3));
  double pyrms = sqrt(abs(sqsum4));

  double ypy = (sumPointStat->_yPy) * den1 - y0*py0;

  //write (z, z0avg * DefConstXl,y0*DefConstXl, yrms*DefConstXl, py0, pyrms, -ypy*DefConstXl, epy*DefConstXl);
}

void write26(ColStat* sumStatZ, ColStat* sumStatPz,  PointStat* sumPointStat, double den1)
{
  double z0  = (sumStatZ->_sum[0]) * den1;
  double sqz = (sumStatZ->_sum[1]) * den1;
  double sqsum5 = sqz - z0*z0;

  double pz0 = (sumStatPz->_sum[0]) * den1;
  double sqpz = (sumStatPz->_sum[1])* den1;
  double sqsum6 = sqpz - pz0*pz0;
  
  double zrms = sqrt(abs(sqsum5));
  double pzrms = sqrt(abs(sqsum6));

  double zpz = (sumPointStat->_zPz) * den1 - z0*pz0;
  //write (z, z0avg * DefConstXl, /*y0*DefConstXl*/ zrms*DefConstXl, pz0, pzrms, -zpz*DefConstXl, epz*DefConstXl);
}

void write27(PointStat* pointStat)
{
  
  //glmax[6] = (pointStat->_max);
  // write (z,z0avg * DefConstXl, glmax(1)*xl,glmax(2),glmax(3)*xl, glmax(4),glmax(5)*xl,glmax(6))

}

void write28()
{
  // npc : max/min particles among cores 
  // nptot: num all particles 
  // write (z,z0avg*xl,npctmin,npctmax,nptot)

}

void write29(ColStat* sumStatX,
	     ColStat* sumStatPx,
	     ColStat* sumStatY, 
	     ColStat* sumStatPy,
	     ColStat* sumStatZ,
	     ColStat* sumStatPz, 
	     double den1)
{
  double x03 = get03(sumStatX, den1);
  double y03 = get03(sumStatY, den1);
  double px03 = get03(sumStatPx, den1);
  double py03 = get03(sumStatPy, den1);


  double z03  = pow((sumStatZ->_sum[2]) * den1,  1.0/3.0);
  double pz03 = pow((sumStatPz->_sum[2]) * den1, 1.0/3.0);
  
  // write(z,z0avg*xl,x03*xl,px03,y03*xl,py03,z03*xl, pz03)
}

void write30(ColStat* sumStatX,
	     ColStat* sumStatPx,
	     ColStat* sumStatY, 
	     ColStat* sumStatPy,
	     ColStat* sumStatZ,
	     ColStat* sumStatPz, 
	     double den1)
{
  double x04 = get04(sumStatX, den1);
  double y04 = get04(sumStatY, den1);
  double px04 = get04(sumStatPx, den1);
  double py04 = get04(sumStatPy, den1);


  double z04  = sqrt(sqrt((sumStatZ->_sum[3]) * den1));
  double pz04 = sqrt(sqrt((sumStatPz->_sum[3]) * den1));

  // write(z,z0avg*xl,x04*xl,px04,y04*xl,py04,z04*xl,pz04)
}

void diag(double z1, double mass, double Scxlt,	  
	  double* data, 
	  int nrows,
	  int Npt) // size of the bunch
{
    double DefConstQmc = mass/1.0e6;
    //double xl = Scxlt; // given
    double DefConstXl = Scxlt;
    double DefConstDen1 = 0.0;

    double xt = 180.0/(2.0*asin(1.0)); // Rad2deg; //given pi = 2.0*asin(1.0) and Rad2deg = 180.0/Pi

    MPI_Comm comm = MPI_COMM_WORLD;
    int i, rank, ierr, j;
    MPI_Comm_rank (comm, &rank);    

    int innp = nrows; 
    int Nbunch = 1;
    uint64_t nptot = 0;  // sum of all pts in all bunchs
    uint64_t innpmb = 0; // sum of local pts accross bunchs in this core
    //z0lc = 0.0;
    //gamlc = 0.0;
    //pz0lc = 0.0;

    double tmp3lc[3], tmp3gl[3];

    int ib=0;
    for (ib=0; i<Nbunch; i++) {
      innpmb = innpmb + innp;
      nptot = nptot + Npt;
      sumtmp3lc(data, nrows, tmp3lc);
    }

	      
    MPI_Allreduce(&tmp3lc, &tmp3gl, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

    DefConstDen1 = 1.0/nptot;
    double z0avg = tmp3gl[0]/nptot;
    double pz0avg = tmp3gl[1]/nptot;
    double gamavg = tmp3gl[2]/nptot;

    /*
    */

    ColStat statX, statY, statZ, statPx, statPy, statPz;
    init(&statX); init(&statY); init(&statZ); 
    init(&statPx); init(&statPy); init(&statPz); 

    PointStat pointStat;
    initPoint(&pointStat);
    for (ib = 1; ib< Nbunch; ib++) {
         for (i = 1; i< innp; i++) {
	      onColumn(data, i, DefColumnX,  &statX,  0.0);   // x0lc, sqsum1local, x0lc3, x0lc4	 
	      onColumn(data, i, DefColumnY,  &statY,  0.0);   // y0lc, sqsum3local, y0lc3, y0lc4 
	      onColumn(data, i, DefColumnPx, &statPx, 0.0);   // px0lc, sqsum2local, px0lc3, px0lc4
	      onColumn(data, i, DefColumnPy, &statPy, 0.0);   // py0lc, sqsum4local, py0lc3, py0lc4

	      // note: z0lc was not really commputed in fortran code
	      onColumn(data, i, DefColumnZ,  &statZ,  z0avg); // z0lc, sqsum5local, z0lc3, z0lc4 
	      onColumn(data, i, DefColumnPz, &statPz, pz0avg);// pz0lc, sqsum6local, pz0lc3, pz0lc4
	      onPoint(data, i,  &pointStat, z0avg, pz0avg, gamavg);
	 }								  
    }

    /*
        tmplc(1) = x0lc;
	tmplc(2) = px0lc;
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc

        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local

        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal

        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4

        tmplc(28) = gamlc
        tmplc(29) = gam2lc
    */  
	

    double glrmax;
    double glmax[6];

    PointStat sumPointStat; initPoint(&sumPointStat);
    SummarizePoint(&pointStat, &sumPointStat);

    ColStat sumX, sumY, sumZ, sumPx, sumPy, sumPz;
    init(&sumX); init(&sumY); init(&sumZ); 
    init(&sumPx); init(&sumPy); init(&sumPz); 

    SummarizeCol(&statX, &sumX); SummarizeCol(&statPx, &sumPx);
    SummarizeCol(&statY, &sumY); SummarizeCol(&statPy, &sumPy);
    SummarizeCol(&statZ, &sumZ); SummarizeCol(&statPz, &sumPz);
    

    //MPI_Reduce(tmplc, tmpgl, 29, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&(pointStat._lcrmax), &glrmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(pointStat._max, glmax, 6, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    double npctmin, npctmax;
    MPI_Reduce(&innpmb, &npctmin, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&innpmb, &npctmax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank == 0) {
      /*
      */
    }
}


