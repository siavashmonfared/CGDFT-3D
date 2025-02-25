#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <stdlib.h>
#include <set>
#include <iterator>
#include <numeric> 
#include <chrono>
#include <random>
#include <complex>
#include <algorithm>
#include <mpi.h>
#include "CGLinearSolver.h"

using namespace std;

#define HALF 13
#define Q    26
#define QLG 6
#define HLG 3

double SD = 1/sqrt(2.0);
double CD = 1/sqrt(3.0);
double KD = 1/sqrt(6.0);
double MD = sqrt(2.0)/sqrt(3.0);
double SL = sqrt(2.0);
double CL = sqrt(3.0);

const double _nx[Q] = {1,0,0,   1, 1,  0,0, 1,1,   1 ,1,1 ,1,   -1,0,0,   -1,-1,0 ,0 ,-1,-1,   -1,-1,-1,-1};
const double _ny[Q] = {0,1,0,  -1, 1,  1,1, 0,0,   1 ,1,-1,-1,  0,-1,0,   1 ,-1,-1,-1,0 ,0 ,   -1,-1,1 ,1 };
const double _nz[Q] = {0,0,1,   0, 0, -1,1,-1,1,   -1,1,1 ,-1,  0,0,-1,   0 ,0 ,1 ,-1,1 ,-1,   1 ,-1,-1,1 };

const double _G2Lxx[Q] = {1., 0., 0.,    0.7071,  0.7071,      0.,      0.,  0.7071, 0.7071, 0.5774, 0.5774, 0.5774, 0.5774,   -1., 0., 0.,   -0.7071,-0.7071,     0.,     0.,-0.70711,-0.7071,-0.5774, -0.5774,-0.5774,-0.5774};
const double _G2Lxy[Q] = {0., 1., 0.,   -0.7071,  0.7071,  0.7071,  0.7071,	     0.,     0., 0.5774, 0.5774,-0.5774,-0.5774,    0.,-1., 0.,    0.7071,-0.7071,-0.7071,-0.7071,	    0.,     0.,-0.5774, -0.5774, 0.5774, 0.5774};
const double _G2Lxz[Q] = {0., 0., 1.,   	 0.,      0., -0.7071,  0.7071, -0.7071, 0.7071,-0.5774, 0.5774, 0.5774,-0.5774,    0., 0.,-1.,         0, 	   0., 0.7071,-0.7071, 0.70711,-0.7071, 0.5774, -0.5774,-0.5774, 0.5774};
const double _G2Lyx[Q] = {0.,-1., 0.,    0.7071, -0.7071,      -1,      -1,	     0.,	 0.,-0.7071,-0.7071, 0.7071, 0.7071,    0., 1., 0.,   -0.7071, 0.7071,      1,      1,	    0.,     0., 0.7071,  0.7071,-0.7071,-0.7071};
const double _G2Lyy[Q] = {1., 0., 1.,    0.7071,  0.7071,      0., 	    0.,	     1.,	 1., 0.7071, 0.7071, 0.7071, 0.7071,   -1., 0., 1.,   -0.7071,-0.7071,      0,      0,	   -1.,	   -1.,-0.7071, -0.7071,-0.7071,-0.7071};
const double _G2Lyz[Q] = {0., 0., 0.,        0.,      0.,      0.,      0.,      0.,	 0.,      0,      0,      0,      0,    0., 0., 0.,        0., 	   0.,      0,	    0,      0.,	    0.,     0.,      0.,     0.,     0.};
const double _G2Lzx[Q] = {0., 0.,-1.,   	 0.,      0.,      0.,      0.,  0.7071,-0.7071,     KD,    -KD,    -KD,     KD,    0., 0., 1.,   	  0., 	   0.,      0,      0, 0.70711,-0.7071,     KD,-0.40825,    -KD,     KD};
const double _G2Lzy[Q] = {0., 0., 0.,   	 0.,      0.,  0.7071, -0.7071,      0.,     0.,     KD,    -KD,	 KD,	-KD,    0., 0., 0.,   	  0., 	   0., 0.7071,-0.7071,      0.,     0.,     KD, -0.4082,     KD,    -KD};
const double _G2Lzz[Q] = {1., 1., 0.,   	 1.,      1.,  0.7071,  0.7071,  0.7071, 0.7071,     MD,     MD,     MD,     MD,    1., 1., 0.,   	  1., 	   1., 0.7071, 0.7071, 0.70711, 0.7071,     MD, 0.81649,     MD,     MD};

const double _L2Gxx[Q] = {1.,  0.,  0.,   0.7071,  0.7071,       0,      0,  0.7071, 0.7071, 0.5774, 0.5774, 0.5774, 0.5774,   -1., 0., 0.,   -0.7071,-0.7071,      0,      0,-0.7071,-0.7071,-0.5774,-0.5774,-0.5774,-0.5774};
const double _L2Gxy[Q] = {0., -1.,  0.,   0.7071, -0.7071,      -1,     -1,	     0.,     0.,-0.7071,-0.7071, 0.7071, 0.7071,    0., 1., 0.,   -0.7071, 0.7071,      1,      1,	    0,     0., 0.7071, 0.7071,-0.7071,-0.7071};
const double _L2Gxz[Q] = {0.,  0., -1.,       0.,      0.,       0,      0,  0.7071,-0.7071, 0.4082,-0.4082,-0.4082, 0.4082,    0., 0., 1.,         0, 	   0.,      0,      0, 0.7071,-0.7071, 0.4082,-0.4082,-0.4082, 0.4082};
const double _L2Gyx[Q] = {0.,  1.,  0.,  -0.7071,  0.7071,  0.7071, 0.7071,	     0.,	 0., 0.5774, 0.5774,-0.5774,-0.5774,    0.,-1., 0.,    0.7071,-0.7071,-0.7071,-0.7071,	    0,	   0.,-0.5774,-0.5774, 0.5774, 0.5774};
const double _L2Gyy[Q] = {1.,  0.,  1.,   0.7071,  0.7071,       0, 	 0,	     1.,	 1., 0.7071, 0.7071, 0.7071, 0.7071,   -1., 0., 1.,   -0.7071,-0.7071,      0,      0,	   -1,	  -1.,-0.7071,-0.7071,-0.7071,-0.7071};
const double _L2Gyz[Q] = {0.,  0.,  0.,       0.,      0.,  0.7071,-0.7071,      0.,	 0., 0.4082,-0.4082, 0.4082,-0.4082,    0., 0., 0.,        0., 	   0., 0.7071,-0.7071,      0,	   0., 0.4082,-0.4082, 0.4082,-0.4082};
const double _L2Gzx[Q] = {0.,  0.,  1.,       0.,      0., -0.7071, 0.7071, -0.7071, 0.7071,-0.5773, 0.5774, 0.5774,-0.5774,    0., 0.,-1.,   	  0., 	   0., 0.7071,-0.7071, 0.7071,-0.7071, 0.5774,-0.5774,-0.5774, 0.5774};
const double _L2Gzy[Q] = {0.,  0.,  0.,       0.,      0.,       0,      0,      0.,	 0.,      0,      0,	  0,      0,    0., 0., 0.,   	  0., 	   0.,      0,      0,      0,	   0.,      0,      0,      0,      0};
const double _L2Gzz[Q] = {1.,  1.,  0.,       1.,      1.,  0.7071, 0.7071,  0.7071, 0.7071, 0.8165, 0.8165, 0.8165, 0.8165,    1., 1., 0.,   	  1., 	   1., 0.7071, 0.7071, 0.7071, 0.7071, 0.8165, 0.8165, 0.8165, 0.8165};

const double _L0[Q] = {1., 1., 1., SL, SL, SL, SL, SL, SL, CL, CL, CL, CL, 1., 1., 1., SL, SL, SL, SL, SL, SL, CL, CL, CL, CL};

struct tensor
{
	double xx, xy, xz;
	double yx, yy, yz;
	double zx, zy, zz;
	void reset() {
		xx = xy = xz = yx = yy = yz = zx = zy = zz = 0.0;
	}
};

struct node { 
	int I[Q]; 
	double IX[Q];
	double IY[Q];
	double IZ[Q];
	size_t ipx[Q];   
	size_t ilg[QLG];
	int nb;  
	int nbo;  
	int nblg;
	int nbolg;
	double FX[Q];
	double FY[Q];
	double FZ[Q];	
	node():nb(0),nbo(0),nblg(0),nbolg(0) { }
};

struct fvec
{
	double x, y, z;
	void reset() {
	x = y = z = 0.0;
	}
};

/*for CG*/
int   ncom;
double *pcom, *xicom, (*nrfunc)(double []);
/*for computation*/
vector<node> nodes;  
vector<size_t> ipImp; 
vector<double> nTr;
vector<size_t> ipTr; 
vector<vector<double> > Ktable; 
double f11, f22, f33;
double* _coef;
double* _coefB;
double* _coefT;
int npore;
/*for MPI*/
MPI_Comm comm3d;
int NeighBor[6];
int S = 0, E = 1, N = 2, W = 3, F = 4, B = 5;
int* blocklen;
int* disp_b;
int* disp_f;
int* blocklen_dt;
int* disp_tp;
int* disp_bt;
int* blocklenID;
int* disp_bID;
int* disp_fID;
int* blocklen_dtID;
int* disp_tpID;
int* disp_btID;
MPI_Datatype tptype, bttype, frtype, bktype, EW_face, BF_face;
MPI_Datatype tptypeID, bttypeID, frtypeID, bktypeID, EW_faceID, BF_faceID;
MPI_Datatype tptypeETA, bttypeETA, frtypeETA, bktypeETA, EW_faceETA, BF_faceETA;
MPI_Status status;
int flag = 1; 
int* lengthX;
int* lengthY;
int* lengthZ;
double* p;
double* force;
int* Id;
int nx_loc, ny_loc, nz_loc, xsum, ysum, zsum, nx_fill, ny_fill, nz_fill, ndof;
/*DFT*/
double* Eta;
double* rho;
double wff,wmf,mu,Temp, Boltz,yvar;
double Lx_lt_shift, Lx_rt_shift, Ly_lt_shift, Ly_rt_shift, Lz_up_shift, Lz_dn_shift;
double* sigF;
/*readInput*/
size_t type, nit;
double Tstar, Tcw;
int nx, ny, nz;	   
int nxp, nyp, nzp;
double gridStep, K1,K0; 	   		   
double strain, ftol, Rpore; 	
double XC, YC, ZC, obj_tol;
double rho_thresh_ub, rho_thresh_lb;
size_t nu_part;
int Lnx, Lny, Lnz;

int ipx(int x, int y, int z, int xsize, int ysize)
{
	return 6 * (x + xsize * y + xsize * ysize * z);
}

void decompose1d (int length[], int ndr ,int npr){

int rem = ndr % npr;

if(rem == 0){
int intL = ndr / npr;
for(int k = 0 ; k < npr; k++) length[k] = intL;
}
	
if (rem == 1){
int intL = (ndr - rem) / npr;
for(int k = 0 ; k < (npr-1); k++) length[k] = intL;
length[npr-1] = rem + intL;
}
	
if (rem > 1){
for(int k = 0 ; k < (npr-1); k++) length[k] = (ndr - rem) / npr;
length[npr-1] = (ndr - rem) / npr + rem ;
}
}


void init(double GridStep, int ndof, size_t nx_loc, size_t ny_loc, size_t nz_loc, int xsum, int ysum, int zsum){

	bool xulF, yulF, xllF, yllF, zulF, zllF;
	bool xulT, yulT, xllT, yllT, zulT, zllT;
	bool xug, xlg, yug, ylg, zug, zlg;
	
	nodes.resize(nx_loc * ny_loc * nz_loc);
	
for (size_t z = 0 ; z < nz_loc ; z++){					
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {
				
				int x_global = xsum + x;
				int y_global = ysum + y;
				int z_global = zsum + z;

				xug = (x_global == nx - 1) ? false : true;
				xlg = (x_global == 0) ? false : true;
				yug = (y_global == ny - 1) ? false : true;
				ylg = (y_global == 0) ? false : true;
				zug = (z_global == nz - 1) ? false : true;
				zlg = (z_global == 0) ? false : true;

				xulF = (x == nx_loc - 1) ? false : true;//if false - continue
				yulF = (y == ny_loc - 1) ? false : true;//if false - continue 
				zulF = (z == nz_loc - 1) ? false : true;//if false - continue
				xllF = (x == 0) ? false : true; //if false - continue 
				yllF = (y == 0) ? false : true; //if false - continue
				zllF = (z == 0) ? false : true; //if false - continue

				xulT = (x == nx_loc - 1) ? true : false;//if true - continue
				yulT = (y == ny_loc - 1) ? true : false;//if true - continue 
				zulT = (z == nz_loc - 1) ? true : false;//if ture - continue
				xllT = (x == 0) ? true : false; //if true - continue 
				yllT = (y == 0) ? true : false; //if true - continue
				zllT = (z == 0) ? true : false; //if true - continue
				
				int iv = 0;
				int ivlg = 0;
				if (xug && xulF){ 
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x+1,y,z,nx_loc,ny_loc);
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[ivlg++] = ipx(x+1,y,z,nx_loc,ny_loc);
				}
				if (xulT){ 
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + z * (ny_loc) + y) * ndof;
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[ivlg++] = (nx_loc * ny_loc * nz_loc + z * (ny_loc) + y) * ndof;
				}
				
				if (yug && yulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x,y+1,z,nx_loc,ny_loc);
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[ivlg++] = ipx(x,y+1,z,nx_loc,ny_loc);
				}
				if (yulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc * 2 + z * (nx_loc + 2) + (x+1)) * ndof;
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[ivlg++] = (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc * 2 + z * (nx_loc + 2) + (x+1)) * ndof;
				}
				
				if (zug && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x,y,z+1,nx_loc,ny_loc);
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[ivlg++] = ipx(x,y,z+1,nx_loc,ny_loc);
				}
				if (zulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * (nx_loc + 2) * nz_loc + (nx_loc + 2) * (ny_loc + 2) + (y + 1) * (nx_loc + 2) + (x + 1)) * ndof;
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[ivlg++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * (nx_loc + 2) * nz_loc + (nx_loc + 2) * (ny_loc + 2) + (y + 1) * (nx_loc + 2) + (x + 1)) * ndof;				
				}
				
				if (xug && ylg){
				if (xulF && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x+1,y-1,z,nx_loc,ny_loc);
				if (yllF && xulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + z * ny_loc + (y-1)) * ndof;
				}
				if (yllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + (nx_loc + 2) * nz_loc + z * (nx_loc + 2) + (x+2)) * ndof;
				}
				}
				
				if (xug && yug){
				if (xulF && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x+1,y+1,z,nx_loc,ny_loc);
				if (xulT && yulF){ 
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + z * (ny_loc) + (y+1) ) * ndof;	
				}
				if (yulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + z * (nx_loc+2) + (x+2)) * ndof;	
				}
				}
				
				if (yug && zlg){
				if (yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x,y+1,z-1,nx_loc,ny_loc);
				if (yulT && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + (z-1) * (nx_loc+2) + (x+1)) * ndof;
				}
				if (zllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * nz_loc * (nx_loc + 2) + (y+2) * (nx_loc+2) + (x+1)) * ndof;
				}
				}
				
				if (yug && zug){
				if (yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x,y+1,z+1,nx_loc,ny_loc);
				if (yulT && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ( nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + (z+1) * (nx_loc + 2) + (x+1)) * ndof;	
				}	
				if (zulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ( nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * nz_loc * (nx_loc + 2) + (nx_loc + 2) * (ny_loc + 2) + (y + 2) * (nx_loc + 2) + (x + 1) ) * ndof;	
				}
				}
				
				if (xug && zlg){
				if (xulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x+1,y,z-1,nx_loc,ny_loc);	
				if (xulT && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + (ny_loc) * (z-1) + y) * ndof;	
				}
				if (zllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * nz_loc * (nx_loc + 2) + (y+1) * (nx_loc + 2) + (x+2) ) * ndof;		
				}
				}
				
				if (xug && zug){
				if (xulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x+1,y,z+1,nx_loc,ny_loc);
				if (xulT && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + (z+1) * ny_loc + y) * ndof;
				}
				if (zulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * nz_loc * (nx_loc + 2) + (nx_loc + 2) * (ny_loc + 2) + (y+1) * (nx_loc + 2) + (2 + x)) * ndof;
				}
				}
				
				if (xug && yug && zlg){
				if (xulF && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x+1,y+1,z-1,nx_loc,ny_loc);
				if (zllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + 2 * nz_loc * (nx_loc + 2) + (y+2) * (nx_loc + 2) + (2 + x)) * ndof;
				}
				if (xulT && yulF && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + (z-1) * ny_loc + (y + 1) ) * ndof;	
				}
				if (yulT && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + (z-1) * (nx_loc + 2) + 1 + (x+1)) * ndof;	
				}
				}
				
				if (xug && yug && zug){
				if (xulF && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x+1,y+1,z+1,nx_loc,ny_loc);
				if (zulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + 2 * nz_loc * (nx_loc + 2) + (ny_loc + 2) * (nx_loc + 2) + (y+2) * (nx_loc + 2) + (2 + x)) * ndof;
				}
				if (xulT && yulF && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + (z+1) * ny_loc + (y+1)) * ndof;	
				}
				if (yulT && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + (z+1) * (nx_loc + 2) + (x + 2)  ) * ndof;	
				}
				}
				
				if (xug && ylg && zug){
				if (xulF && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x+1,y-1,z+1,nx_loc,ny_loc);
				if (zulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + 2 * nz_loc * (nx_loc + 2) + (ny_loc + 2) * (nx_loc + 2) + (y) * (nx_loc + 2) + (2 + x)) * ndof;
				}
				if (xulT && yllF && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + (z+1) * ny_loc + (y - 1)  ) * ndof;	
				}
				if (yllT && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + nz_loc * (nx_loc + 2) + (z+1) * (nx_loc + 2) + (x + 2)  ) * ndof;
				}
				}
				
				if (xug && ylg && zlg){
				if (xulF && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = ipx(x+1,y-1,z-1,nx_loc,ny_loc);
				if (zllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + 2 * nz_loc * (nx_loc + 2) + y * (nx_loc + 2) + (x + 2)  ) * ndof;	
				}
				if (xulT && yllF && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + (z-1) * ny_loc + (y - 1) ) * ndof;	
				}
				if (yllT && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[iv++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + nz_loc * (nx_loc + 2) + (z-1) * (nx_loc + 2) + (x + 2)    ) * ndof;	
				}
				}
	
				int ivo = 0;
				int ivolg = 0;
				if (xlg && xllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = ipx(x-1,y,z,nx_loc,ny_loc);
				if (xlg && xllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[HLG+ivolg++] = ipx(x-1,y,z,nx_loc,ny_loc);
				if (xllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc + z * ny_loc + (y)) * ndof;
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[HLG+ivolg++] = (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc + z * ny_loc + (y)) * ndof;
				}

				if (ylg && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = ipx(x, y - 1,z,nx_loc,ny_loc);
				if (ylg && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[HLG+ivolg++] = ipx(x, y - 1,z,nx_loc,ny_loc);
				if (yllT){ 
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc * 2 + (nx_loc + 2) * nz_loc + z * (nx_loc + 2) + (x+1)) * ndof;
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[HLG+ivolg++] = (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc * 2 + (nx_loc + 2) * nz_loc + z * (nx_loc + 2) + (x+1)) * ndof;
				}
				
				if (zlg && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = ipx(x,y,z-1,nx_loc,ny_loc);
				if (zlg && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[HLG+ivolg++] = ipx(x,y,z-1,nx_loc,ny_loc);
				if (zllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * (nx_loc + 2) * nz_loc + (y+1) * (nx_loc + 2) + (x+1)) * ndof;
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ilg[HLG+ivolg++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * (nx_loc + 2) * nz_loc + (y+1) * (nx_loc + 2) + (x+1)) * ndof;
				}
					
				if (xlg && yug){
				if (xllF && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = ipx(x-1,y+1,z,nx_loc,ny_loc);
				if (xllT && yulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + nz_loc * ny_loc + z * (ny_loc) + (y+1)) * ndof;
				}
				if (yulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + z * (nx_loc + 2) + (x)) * ndof;
				}
				}
					
				if (xlg && ylg){
				if (xllF && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = ipx(x-1,y-1,z,nx_loc,ny_loc);
				if (xllT && yllF){ 
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc + z * (ny_loc) + (y-1)) * ndof;
				}
				if (yllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + (nx_loc + 2) * nz_loc + z * (nx_loc + 2) + (x)) * ndof;	
				}
				}
				
				if (ylg && zug){
				if (yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = ipx(x,y-1,z+1,nx_loc,ny_loc);
				if (yllT && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + nz_loc * (nx_loc + 2) + (z+1) * (nx_loc + 2) + (1 + x)) * ndof;
				}
				if (zulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + 2 * nz_loc * (nx_loc + 2) + (nx_loc + 2) * (ny_loc + 2) + (y) * (nx_loc + 2) + (1 + x) ) * ndof;	
				}
				}
				
				if (ylg && zlg){
				if (yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = ipx(x,y-1,z-1,nx_loc,ny_loc);
				if (yllT && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + nz_loc * (nx_loc + 2) + (z-1) * (nx_loc + 2) + 1 + (x)) * ndof;
				}
				if (zllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * nz_loc * (nx_loc + 2) + y * (nx_loc + 2) + 1 + (x)) * ndof;
				}
				}
				
				if (xlg && zug){
				if (xllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = ipx(x-1,y,z+1,nx_loc,ny_loc);
				if (xllT && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc + (z+1) * ny_loc + (y)) * ndof;	
				}
				if (zulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * nz_loc * (nx_loc + 2) + (nx_loc + 2) * (ny_loc + 2) + (y+1) * (nx_loc + 2) + (x)) * ndof;	
				}
				}
				
				if (xlg && zlg){
				if (xllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = ipx(x-1,y,z-1,nx_loc,ny_loc);
				if (xllT && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + nz_loc * ny_loc + (z-1) * ny_loc + (y)) * ndof;	
				}
				if (zllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF+ivo++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + 2 * nz_loc * (nx_loc + 2) + (y+1) * (nx_loc + 2) + (x)) * ndof;	
				}
				}
								
				if (xlg && ylg && zug){
				if (xllF && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = ipx(x-1,y-1,z+1,nx_loc,ny_loc);
				if (zulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + 2 * nz_loc * (nx_loc + 2) + (nx_loc + 2) * (ny_loc + 2) + (y) * (nx_loc + 2) + (x)) * ndof;	
				}
				if (xllT && yllF && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc + (z+1) * ny_loc + (y - 1) ) * ndof;		
				}
				if (yllT && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + nz_loc * (nx_loc + 2) + (z+1) * (nx_loc + 2) + (x)) * ndof;		
				}
				}
								
				if (xlg && ylg && zlg){
				if (xllF && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = ipx(x-1,y-1,z-1,nx_loc,ny_loc);
				if (zllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + 2 * nz_loc * (nx_loc + 2) + (y) * (nx_loc + 2) + (x)) * ndof;	
				}
				if (xllT && yllF && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc + (z-1) * ny_loc + (y - 1) ) * ndof;		
				}
				if (yllT && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + nz_loc * (nx_loc + 2) + (z-1) * (nx_loc + 2) + (x)) * ndof;		
				}
				}
				
				if (xlg && yug && zlg){
				if (xllF && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = ipx(x-1,y+1,z-1,nx_loc,ny_loc);
				if (zllT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + 2 * nz_loc * (nx_loc + 2) + (y + 2) * (nx_loc + 2) + (x) ) * ndof;
				}	
				if (xllT && yulF && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + nz_loc * ny_loc + (z-1) * ny_loc + (y + 1) ) * ndof;
				}
				if (yulT && zllF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + (z-1) * (nx_loc + 2) + (x) ) * ndof;
				}
				}
				
				if (xlg && yug && zug){
				if (xllF && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = ipx(x-1,y+1,z+1,nx_loc,ny_loc);
				if (zulT){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + 2 * nz_loc * (nx_loc + 2) + (nx_loc + 2) * (ny_loc + 2) + (y + 2) * (nx_loc + 2) + (x) ) * ndof;
				}	
				if (xllT && yulF && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + nz_loc * ny_loc + (z+1) * ny_loc + (y + 1) ) * ndof;
				}
				if (yulT && zulF){
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[HALF + ivo++] = (nx_loc * ny_loc * nz_loc + 2 * nz_loc * ny_loc + (z+1) * (nx_loc + 2) + (x) ) * ndof;
				}
				}

				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nb = iv;
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nbo = ivo;
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nblg = ivlg;
				nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nbolg = ivolg;
				
				
				iv = 0;
				if (xug && xulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 0;
				if (xug && xulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 0;
				
				if (yug && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 1;
				if (yug && yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 1;
				
				if (zug && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 2;
				if (zug && zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 2;
				
				if (xug && ylg){
				if (xulF && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 3;
				if (yllF && xulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 3;
				if (yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 3;
				}
				
				if (xug && yug){
				if (xulF && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 4;
				if (xulT && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 4;
				if (yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 4;
				}
				
				if (yug && zlg){
				if (yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 5;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 5;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 5;
				}
				
				if (yug && zug){
				if (yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 6;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 6;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 6;
				}
				
				if (xug && zlg){
				if (xulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 7;	
				if (xulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 7;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 7;
				}
				
				if (xug && zug){
				if (xulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 8;
				if (xulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 8;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 8;
				}
				
				if (xug && yug && zlg){
				if (xulF && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 9;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 9;
				if (xulT && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 9;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 9;
				}
				
				if (xug && yug && zug){
				if (xulF && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 10;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 10;
				if (xulT && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 10;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 10;
				}
				
				if (xug && ylg && zug){
				if (xulF && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 11;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 11;
				if (xulT && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 11;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 11;
				}
				
				if (xug && ylg && zlg){
				if (xulF && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 12;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 12;
				if (xulT && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 12;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[iv++] = 12;
				}
	
				ivo = 0;
				if (xlg && xllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF + ivo++] = 13;
				if (xlg && xllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF + ivo++] = 13;

				if (ylg && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF + ivo++] = 14;
				if (ylg && yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF + ivo++] = 14;
				
				if (zlg && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF + ivo++] = 15;
				if (zlg && zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF + ivo++] = 15;
					
				if (xlg && yug){
				if (xllF && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 16;
				if (xllT && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 16;
				if (yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 16;
				}
					
				if (xlg && ylg){
				if (xllF && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 17;
				if (xllT && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 17;
				if (yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 17;
				}
				
				if (ylg && zug){
				if (yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 18;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 18;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 18;
				}
				
				if (ylg && zlg){
				if (yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 19;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 19;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 19;
				}
				
				if (xlg && zug){
				if (xllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 20;
				if (xllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 20;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 20;
				}
				
				if (xlg && zlg){
				if (xllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 21;	
				if (xllT && zllF)  nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 21;	
				if (zllT)  nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 21;	
				}
								
				if (xlg && ylg && zug){
				if (xllF && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 22;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 22;
				if (xllT && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 22;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 22;
				}
								
				if (xlg && ylg && zlg){
				if (xllF && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 23;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 23;
				if (xllT && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 23;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 23;
				}
				
				if (xlg && yug && zlg){
				if (xllF && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 24;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 24;
				if (xllT && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 24;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 24;
				}
				
				if (xlg && yug && zug){
				if (xllF && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 25;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 25;
				if (xllT && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 25;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[HALF+ivo++] = 25;
				}
				
				iv = 0;
				if (xug && xulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (xug && xulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				
				if (yug && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = x_global * gridStep;
				if (yug && yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = x_global * gridStep;
				
				if (zug && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = x_global * gridStep;
				if (zug && zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = x_global * gridStep;
				
				if (xug && ylg){
				if (xulF && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (yllF && xulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				}
				
				if (xug && yug){
				if (xulF && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (xulT && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				}
				
				if (yug && zlg){
				if (yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = x_global * gridStep;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = x_global * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = x_global * gridStep;
				}
				
				if (yug && zug){
				if (yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = x_global * gridStep;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = x_global * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = x_global * gridStep;
				}
				
				if (xug && zlg){
				if (xulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;	
				if (xulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				}
				
				if (xug && zug){
				if (xulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (xulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				}
				
				if (xug && yug && zlg){
				if (xulF && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (xulT && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				}
				
				if (xug && yug && zug){
				if (xulF && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (xulT && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				}
				
				if (xug && ylg && zug){
				if (xulF && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (xulT && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				}
				
				if (xug && ylg && zlg){
				if (xulF && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (xulT && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[iv++] = (x_global+1) * gridStep;
				}
	
				ivo = 0;
				if (xlg && xllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF + ivo++] = (x_global-1) * gridStep;
				if (xlg && xllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF + ivo++] = (x_global-1) * gridStep;

				if (ylg && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF + ivo++] = x_global * gridStep;
				if (ylg && yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF + ivo++] = x_global * gridStep;
				
				if (zlg && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF + ivo++] = x_global * gridStep;
				if (zlg && zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF + ivo++] = x_global * gridStep;
					
				if (xlg && yug){
				if (xllF && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (xllT && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				}
					
				if (xlg && ylg){
				if (xllF && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (xllT && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				}
				
				if (ylg && zug){
				if (yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = x_global * gridStep;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = x_global * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = x_global * gridStep;
				}
				
				if (ylg && zlg){
				if (yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = x_global * gridStep;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = x_global * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = x_global * gridStep;
				}
				
				if (xlg && zug){
				if (xllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (xllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				}
				
				if (xlg && zlg){
				if (xllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;	
				if (xllT && zllF)  nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;	
				if (zllT)  nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;	
				}
								
				if (xlg && ylg && zug){
				if (xllF && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (xllT && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				}
								
				if (xlg && ylg && zlg){
				if (xllF && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (xllT && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				}
				
				if (xlg && yug && zlg){
				if (xllF && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (xllT && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				}
				
				if (xlg && yug && zug){
				if (xllF && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (xllT && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IX[HALF+ivo++] = (x_global-1) * gridStep;
				}
				
				iv = 0;
				if (xug && xulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global) * gridStep;
				if (xug && xulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global) * gridStep;
				
				if (yug && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				if (yug && yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				
				if (zug && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global) * gridStep;
				if (zug && zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global) * gridStep;
				
				if (xug && ylg){
				if (xulF && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global-1) * gridStep;
				if (yllF && xulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global-1) * gridStep;
				if (yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global-1) * gridStep;
				}
				
				if (xug && yug){
				if (xulF && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				if (xulT && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				if (yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				}
				
				if (yug && zlg){
				if (yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				}
				
				if (yug && zug){
				if (yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				}
				
				if (xug && zlg){
				if (xulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global) * gridStep;	
				if (xulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global) * gridStep;
				}
				
				if (xug && zug){
				if (xulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global) * gridStep;
				if (xulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global) * gridStep;
				}
				
				if (xug && yug && zlg){
				if (xulF && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global + 1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global + 1) * gridStep;
				if (xulT && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global + 1) * gridStep;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global + 1) * gridStep;
				}
				
				if (xug && yug && zug){
				if (xulF && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				if (xulT && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global+1) * gridStep;
				}
				
				if (xug && ylg && zug){
				if (xulF && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global-1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global-1) * gridStep;
				if (xulT && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global-1) * gridStep;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global-1) * gridStep;
				}
				
				if (xug && ylg && zlg){
				if (xulF && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global-1) * gridStep;
				if (xulT && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global-1) * gridStep;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[iv++] = (y_global-1) * gridStep;
				}
	
				ivo = 0;
				if (xlg && xllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF + ivo++] = (y_global) * gridStep;
				if (xlg && xllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF + ivo++] = (y_global) * gridStep;

				if (ylg && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF + ivo++] = (y_global-1) * gridStep;
				if (ylg && yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF + ivo++] = (y_global-1) * gridStep;
				
				if (zlg && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF + ivo++] = y_global * gridStep;
				if (zlg && zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF + ivo++] = y_global * gridStep;
					
				if (xlg && yug){
				if (xllF && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global+1) * gridStep;
				if (xllT && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global+1) * gridStep;
				if (yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global+1) * gridStep;
				}
					
				if (xlg && ylg){
				if (xllF && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (xllT && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				}
				
				if (ylg && zug){
				if (yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				}
				
				if (ylg && zlg){
				if (yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				}
				
				if (xlg && zug){
				if (xllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global) * gridStep;
				if (xllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global) * gridStep;
				}
				
				if (xlg && zlg){
				if (xllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global) * gridStep;	
				if (xllT && zllF)  nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global) * gridStep;	
				if (zllT)  nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global) * gridStep;	
				}
								
				if (xlg && ylg && zug){
				if (xllF && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (xllT && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				}
								
				if (xlg && ylg && zlg){
				if (xllF && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (xllT && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global-1) * gridStep;
				}
				
				if (xlg && yug && zlg){
				if (xllF && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global+1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global+1) * gridStep;
				if (xllT && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global+1) * gridStep;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global+1) * gridStep;
				}
				
				if (xlg && yug && zug){
				if (xllF && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global+1) * gridStep;
				if (xllT && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global+1) * gridStep;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IY[HALF+ivo++] = (y_global+1) * gridStep;
				}
				
				iv = 0;
				if (xug && xulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global) * gridStep;
				if (xug && xulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global) * gridStep;
				
				if (yug && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global) * gridStep;
				if (yug && yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global) * gridStep;
				
				if (zug && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				if (zug && zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				
				if (xug && ylg){
				if (xulF && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global) * gridStep;
				if (yllF && xulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global) * gridStep;
				if (yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global) * gridStep;
				}
				
				if (xug && yug){
				if (xulF && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global) * gridStep;
				if (xulT && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global) * gridStep;
				if (yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global) * gridStep;
				}
				
				if (yug && zlg){
				if (yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				}
				
				if (yug && zug){
				if (yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				}
				
				if (xug && zlg){
				if (xulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;	
				if (xulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				}
				
				if (xug && zug){
				if (xulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				if (xulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				}
				
				if (xug && yug && zlg){
				if (xulF && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				if (xulT && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				}
				
				if (xug && yug && zug){
				if (xulF && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				if (xulT && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				}
				
				if (xug && ylg && zug){
				if (xulF && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				if (xulT && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global+1) * gridStep;
				}
				
				if (xug && ylg && zlg){
				if (xulF && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				if (xulT && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[iv++] = (z_global-1) * gridStep;
				}
	
				ivo = 0;
				if (xlg && xllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF + ivo++] = (z_global) * gridStep;
				if (xlg && xllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF + ivo++] = (z_global) * gridStep;

				if (ylg && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF + ivo++] = (z_global) * gridStep;
				if (ylg && yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF + ivo++] = (z_global) * gridStep;
				
				if (zlg && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF + ivo++] = (z_global-1) * gridStep;
				if (zlg && zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF + ivo++] = (z_global-1) * gridStep;
					
				if (xlg && yug){
				if (xllF && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global) * gridStep;
				if (xllT && yulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global) * gridStep;
				if (yulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global) * gridStep;
				}
					
				if (xlg && ylg){
				if (xllF && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global) * gridStep;
				if (xllT && yllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global) * gridStep;
				if (yllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global) * gridStep;
				}
				
				if (ylg && zug){
				if (yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				}
				
				if (ylg && zlg){
				if (yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;
				}
				
				if (xlg && zug){
				if (xllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				if (xllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				}
				
				if (xlg && zlg){
				if (xllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;	
				if (xllT && zllF)  nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;	
				if (zllT)  nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;	
				}
								
				if (xlg && ylg && zug){
				if (xllF && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				if (xllT && yllF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				if (yllT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				}
								
				if (xlg && ylg && zlg){
				if (xllF && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;
				if (xllT && yllF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;
				if (yllT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;
				}
				
				if (xlg && yug && zlg){
				if (xllF && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;
				if (zllT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;
				if (xllT && yulF && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;
				if (yulT && zllF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global-1) * gridStep;
				}
				
				if (xlg && yug && zug){
				if (xllF && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				if (zulT) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				if (xllT && yulF && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				if (yulT && zulF) nodes[ipx(x,y,z,nx_loc,ny_loc)/6].IZ[HALF+ivo++] = (z_global+1) * gridStep;
				}
				
			}
		}
	}
}

void muVT_bc(size_t nx_loc, size_t ny_loc, size_t nz_loc, int xsum, int ysum, int zsum){
	
	for (size_t z = 0 ; z < nz_loc ; z++){					
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {
				
				int x_global = xsum + x;
				int y_global = ysum + y;
				int z_global = zsum + z;
				
				if (x_global == 0){
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc));
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 1);
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 2);
				}
				if (x_global == nx - 1){
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc));
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 1);					
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 2);					
				}
				if (y_global == 0){
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc));
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 1);
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 2);
				}
				if (y_global == ny - 1){
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc));
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 1);					
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 2);					
				}
				if (z_global == 0){
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc));
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 1);
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 2);
				}
				if (z_global == nz - 1){
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc));
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 1);					
				ipImp.push_back(ipx(x,y,z,nx_loc,ny_loc) + 2);					
				}
				
				
			}
		}
	}		
}

void readSimParam()
{
    	const char * name = "input_12102018.txt";
    	ifstream file(name);
    	file >> Lnx >> Lny >> Lnz;
	file >> nxp >> nyp >> nzp;
    	file >> XC >> YC >> ZC >> Rpore;
	file >> type >> strain >> nit >> ftol;
	file >> Tstar >> gridStep >> Tcw >> yvar;
	file >> obj_tol >> nu_part;
	file >> rho_thresh_ub >> rho_thresh_lb;
	file >> Lx_lt_shift >> Lx_rt_shift;
	file >> Ly_lt_shift >> Ly_rt_shift;
	file >> Lz_dn_shift >> Lz_up_shift;
}

double Up_springs(double p[]){

	double Up = 0.0;

	size_t ip, jp;
	double kn,kb,kt;
	double dxi,dxj,dyi,dyj,dzi,dzj;
	double g2lxx,g2lxy,g2lxz,g2lyx,g2lyy,g2lyz,g2lzx,g2lzy,g2lzz;
	double tzi,tzj,tyi,tyj,txi,txj;
	double dlx, dly, dlz;
	int p1, p2;

	for (size_t z = 0 ; z < nz_loc ; z++) {
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {
				
			int x_global = xsum + x;
			int y_global = ysum + y;
			int z_global = zsum + z;
				
			size_t k = ipx(x,y,z,nx_loc,ny_loc)/6;
			ip = ipx(x,y,z,nx_loc,ny_loc);
			
			for (int l = 0 ; l < nodes[k].nb ; ++l) {
			jp = nodes[k].ipx[l];

			dxi = p[ip] - 	x_global * gridStep;
			dxj = p[jp] - 	nodes[k].IX[l];
			dyi = p[ip+1] - y_global * gridStep;
			dyj = p[jp+1] - nodes[k].IY[l];
			dzi = p[ip+2] - z_global * gridStep;
			dzj = p[jp+2] - nodes[k].IZ[l];

			txi = p[ip + 3];
			txj = p[jp + 3];
			tyi = p[ip + 4];
			tyj = p[jp + 4];
			tzi = p[ip + 5];
			tzj = p[jp + 5];
			
			kn = _coef[nodes[k].I[l]] * Ktable[Id[k]][Id[jp/6]];
			kb = _coefB[nodes[k].I[l]] * Ktable[Id[k]][Id[jp/6]];
			kt = _coefT[nodes[k].I[l]] * Ktable[Id[k]][Id[jp/6]];
			
			g2lxx = _G2Lxx[nodes[k].I[l]];
			g2lxy = _G2Lxy[nodes[k].I[l]];
			g2lxz = _G2Lxz[nodes[k].I[l]];
			g2lyx = _G2Lyx[nodes[k].I[l]];
			g2lyy = _G2Lyy[nodes[k].I[l]];
			g2lyz = _G2Lyz[nodes[k].I[l]];
			g2lzx = _G2Lzx[nodes[k].I[l]];
			g2lzy = _G2Lzy[nodes[k].I[l]];
			g2lzz = _G2Lzz[nodes[k].I[l]];			
					
			Up += 	(0.1667) * (3*pow(dyi,2)*pow(g2lxy,2)*kn - 6*dyi*dyj*pow(g2lxy,2)*kn + 3*pow(dyj,2)*pow(g2lxy,2)*kn + 
					6*dyi*dzi*g2lxy*g2lxz*kn - 6*dyj*dzi*g2lxy*g2lxz*kn - 
					6*dyi*dzj*g2lxy*g2lxz*kn + 6*dyj*dzj*g2lxy*g2lxz*kn + 
					3*pow(dzi,2)*pow(g2lxz,2)*kn - 6*dzi*dzj*pow(g2lxz,2)*kn + 3*pow(dzj,2)*pow(g2lxz,2)*kn + 
					3*pow(dyi,2)*pow(g2lyy,2)*kt - 6*dyi*dyj*pow(g2lyy,2)*kt + 3*pow(dyj,2)*pow(g2lyy,2)*kt + 
					6*dyi*dzi*g2lyy*g2lyz*kt - 6*dyj*dzi*g2lyy*g2lyz*kt - 
					6*dyi*dzj*g2lyy*g2lyz*kt + 6*dyj*dzj*g2lyy*g2lyz*kt + 
					3*pow(dzi,2)*pow(g2lyz,2)*kt - 6*dzi*dzj*pow(g2lyz,2)*kt + 3*pow(dzj,2)*pow(g2lyz,2)*kt + 
					3*pow(dyi,2)*pow(g2lzy,2)*kt - 6*dyi*dyj*pow(g2lzy,2)*kt + 3*pow(dyj,2)*pow(g2lzy,2)*kt + 
					6*dyi*dzi*g2lzy*g2lzz*kt - 6*dyj*dzi*g2lzy*g2lzz*kt - 
					6*dyi*dzj*g2lzy*g2lzz*kt + 6*dyj*dzj*g2lzy*g2lzz*kt + 
					3*pow(dzi,2)*pow(g2lzz,2)*kt - 6*dzi*dzj*pow(g2lzz,2)*kt + 3*pow(dzj,2)*pow(g2lzz,2)*kt + 
					3*pow(dxi,2)*(pow(g2lxx,2)*kn + (pow(g2lyx,2) + pow(g2lzx,2))*kt) + 
					3*pow(dxj,2)*(pow(g2lxx,2)*kn + (pow(g2lyx,2) + pow(g2lzx,2))*kt) + 
					3*dyi*g2lyy*g2lzx*kt*_L0[nodes[k].I[l]]*txi - 3*dyj*g2lyy*g2lzx*kt*_L0[nodes[k].I[l]]*txi + 
					3*dzi*g2lyz*g2lzx*kt*_L0[nodes[k].I[l]]*txi - 3*dzj*g2lyz*g2lzx*kt*_L0[nodes[k].I[l]]*txi - 
					3*dyi*g2lyx*g2lzy*kt*_L0[nodes[k].I[l]]*txi + 3*dyj*g2lyx*g2lzy*kt*_L0[nodes[k].I[l]]*txi - 
					3*dzi*g2lyx*g2lzz*kt*_L0[nodes[k].I[l]]*txi + 3*dzj*g2lyx*g2lzz*kt*_L0[nodes[k].I[l]]*txi + 
					pow(g2lyx,2)*kt*pow(_L0[nodes[k].I[l]],2)*pow(txi,2) + pow(g2lzx,2)*kt*pow(_L0[nodes[k].I[l]],2)*pow(txi,2) + 
					3*dyi*g2lyy*g2lzx*kt*_L0[nodes[k].I[l]]*txj - 3*dyj*g2lyy*g2lzx*kt*_L0[nodes[k].I[l]]*txj + 
					3*dzi*g2lyz*g2lzx*kt*_L0[nodes[k].I[l]]*txj - 3*dzj*g2lyz*g2lzx*kt*_L0[nodes[k].I[l]]*txj - 
					3*dyi*g2lyx*g2lzy*kt*_L0[nodes[k].I[l]]*txj + 3*dyj*g2lyx*g2lzy*kt*_L0[nodes[k].I[l]]*txj - 
					3*dzi*g2lyx*g2lzz*kt*_L0[nodes[k].I[l]]*txj + 3*dzj*g2lyx*g2lzz*kt*_L0[nodes[k].I[l]]*txj + 
					pow(g2lyx,2)*kt*pow(_L0[nodes[k].I[l]],2)*txi*txj + pow(g2lzx,2)*kt*pow(_L0[nodes[k].I[l]],2)*txi*txj + 
					pow(g2lyx,2)*kt*pow(_L0[nodes[k].I[l]],2)*pow(txj,2) + pow(g2lzx,2)*kt*pow(_L0[nodes[k].I[l]],2)*pow(txj,2) + 
					3*dzi*g2lyz*g2lzy*kt*_L0[nodes[k].I[l]]*tyi - 3*dzj*g2lyz*g2lzy*kt*_L0[nodes[k].I[l]]*tyi - 
					3*dzi*g2lyy*g2lzz*kt*_L0[nodes[k].I[l]]*tyi + 3*dzj*g2lyy*g2lzz*kt*_L0[nodes[k].I[l]]*tyi + 
					2*g2lyx*g2lyy*kt*pow(_L0[nodes[k].I[l]],2)*txi*tyi + 2*g2lzx*g2lzy*kt*pow(_L0[nodes[k].I[l]],2)*txi*tyi + 
					g2lyx*g2lyy*kt*pow(_L0[nodes[k].I[l]],2)*txj*tyi + g2lzx*g2lzy*kt*pow(_L0[nodes[k].I[l]],2)*txj*tyi + 
					pow(g2lyy,2)*kt*pow(_L0[nodes[k].I[l]],2)*pow(tyi,2) + pow(g2lzy,2)*kt*pow(_L0[nodes[k].I[l]],2)*pow(tyi,2) + 
					3*dzi*g2lyz*g2lzy*kt*_L0[nodes[k].I[l]]*tyj - 3*dzj*g2lyz*g2lzy*kt*_L0[nodes[k].I[l]]*tyj - 
					3*dzi*g2lyy*g2lzz*kt*_L0[nodes[k].I[l]]*tyj + 3*dzj*g2lyy*g2lzz*kt*_L0[nodes[k].I[l]]*tyj + 
					g2lyx*g2lyy*kt*pow(_L0[nodes[k].I[l]],2)*txi*tyj + g2lzx*g2lzy*kt*pow(_L0[nodes[k].I[l]],2)*txi*tyj + 
					2*g2lyx*g2lyy*kt*pow(_L0[nodes[k].I[l]],2)*txj*tyj + 2*g2lzx*g2lzy*kt*pow(_L0[nodes[k].I[l]],2)*txj*tyj + 
					pow(g2lyy,2)*kt*pow(_L0[nodes[k].I[l]],2)*tyi*tyj + pow(g2lzy,2)*kt*pow(_L0[nodes[k].I[l]],2)*tyi*tyj + 
					pow(g2lyy,2)*kt*pow(_L0[nodes[k].I[l]],2)*pow(tyj,2) + pow(g2lzy,2)*kt*pow(_L0[nodes[k].I[l]],2)*pow(tyj,2) - 
					3*dyi*g2lyz*g2lzy*kt*_L0[nodes[k].I[l]]*tzi + 3*dyj*g2lyz*g2lzy*kt*_L0[nodes[k].I[l]]*tzi + 
					3*dyi*g2lyy*g2lzz*kt*_L0[nodes[k].I[l]]*tzi - 3*dyj*g2lyy*g2lzz*kt*_L0[nodes[k].I[l]]*tzi + 
					2*g2lyx*g2lyz*kt*pow(_L0[nodes[k].I[l]],2)*txi*tzi + 2*g2lzx*g2lzz*kt*pow(_L0[nodes[k].I[l]],2)*txi*tzi + 
					g2lyx*g2lyz*kt*pow(_L0[nodes[k].I[l]],2)*txj*tzi + g2lzx*g2lzz*kt*pow(_L0[nodes[k].I[l]],2)*txj*tzi + 
					2*g2lyy*g2lyz*kt*pow(_L0[nodes[k].I[l]],2)*tyi*tzi + 2*g2lzy*g2lzz*kt*pow(_L0[nodes[k].I[l]],2)*tyi*tzi + 
					g2lyy*g2lyz*kt*pow(_L0[nodes[k].I[l]],2)*tyj*tzi + g2lzy*g2lzz*kt*pow(_L0[nodes[k].I[l]],2)*tyj*tzi + 
					pow(g2lyz,2)*kt*pow(_L0[nodes[k].I[l]],2)*pow(tzi,2) + pow(g2lzz,2)*kt*pow(_L0[nodes[k].I[l]],2)*pow(tzi,2) + 
					kt*_L0[nodes[k].I[l]]*(3*dyj*(g2lyz*g2lzy - g2lyy*g2lzz) + 
					dyi*(-3*g2lyz*g2lzy + 3*g2lyy*g2lzz) + 
					_L0[nodes[k].I[l]]*(g2lyx*g2lyz*(txi + 2*txj) + 
					g2lzx*g2lzz*(txi + 2*txj) + (g2lyy*g2lyz + 
					g2lzy*g2lzz)*(tyi + 2*tyj) + (pow(g2lyz,2) + 
					pow(g2lzz,2))*tzi))*tzj + (pow(g2lyz,2) + pow(g2lzz,2))*kt*pow(_L0[nodes[k].I[l]],2)*pow(tzj,2) - 
					3*dxj*(2*dzi*g2lxx*g2lxz*kn - 2*dzj*g2lxx*g2lxz*kn + 
					2*dzi*g2lyx*g2lyz*kt - 2*dzj*g2lyx*g2lyz*kt + 
					2*dzi*g2lzx*g2lzz*kt - 2*dzj*g2lzx*g2lzz*kt + 
					2*dyi*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) - 
					2*dyj*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) - 
					g2lyy*g2lzx*kt*_L0[nodes[k].I[l]]*tyi + g2lyx*g2lzy*kt*_L0[nodes[k].I[l]]*tyi - 
					g2lyy*g2lzx*kt*_L0[nodes[k].I[l]]*tyj + g2lyx*g2lzy*kt*_L0[nodes[k].I[l]]*tyj - 
					g2lyz*g2lzx*kt*_L0[nodes[k].I[l]]*tzi + g2lyx*g2lzz*kt*_L0[nodes[k].I[l]]*tzi - 
					g2lyz*g2lzx*kt*_L0[nodes[k].I[l]]*tzj + g2lyx*g2lzz*kt*_L0[nodes[k].I[l]]*tzj) - 
					3*dxi*(2*dyj*g2lxx*g2lxy*kn - 2*dzi*g2lxx*g2lxz*kn + 
					2*dzj*g2lxx*g2lxz*kn + 2*dyj*g2lyx*g2lyy*kt - 
					2*dzi*g2lyx*g2lyz*kt + 2*dzj*g2lyx*g2lyz*kt + 
					2*dyj*g2lzx*g2lzy*kt - 2*dzi*g2lzx*g2lzz*kt + 
					2*dzj*g2lzx*g2lzz*kt + 
					2*dxj*(pow(g2lxx,2)*kn + (pow(g2lyx,2) + pow(g2lzx,2))*kt) - 
					2*dyi*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) + 
					g2lyy*g2lzx*kt*_L0[nodes[k].I[l]]*tyi - g2lyx*g2lzy*kt*_L0[nodes[k].I[l]]*tyi + 
					g2lyy*g2lzx*kt*_L0[nodes[k].I[l]]*tyj - g2lyx*g2lzy*kt*_L0[nodes[k].I[l]]*tyj + 
					g2lyz*g2lzx*kt*_L0[nodes[k].I[l]]*tzi - 
					g2lyx*g2lzz*kt*_L0[nodes[k].I[l]]*tzi + (g2lyz*g2lzx - g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*tzj));
					
					
				}
			}
		}
	}
			
		for (size_t z = 0 ; z < nz_loc ; z++) {					
				for (size_t y = 0 ; y < ny_loc ; y++) {
					for (size_t x = 0 ; x < nx_loc ; x++) {
			
			int x_global = xsum + x;
			int y_global = ysum + y;
			int z_global = zsum + z;
				
			dlx = p[ipx(x,y,z,nx_loc,ny_loc)] - x_global * gridStep;
			dly = p[ipx(x,y,z,nx_loc,ny_loc) + 1] - y_global * gridStep;
			dlz = p[ipx(x,y,z,nx_loc,ny_loc) + 2] - z_global * gridStep;
				
			Up += force[ipx(x,y,z,nx_loc,ny_loc)] * dlx;
			Up += force[ipx(x,y,z,nx_loc,ny_loc) + 1] * dly;
			Up += force[ipx(x,y,z,nx_loc,ny_loc) + 2] * dlz;

					}
				}
			}

	return Up;
}


void Grad_Up_springs_app (double p[], double df[]){

	for (size_t i = 0 ; i < (nx_fill * ny_fill * nz_fill * ndof) ; i++) df[i] = 0.0;

	size_t ip, jpx, itx, jtx;
	size_t ipy, jpy, ity, jty;
	size_t ipz, jpz, itz, jtz;
	
	double dxi,dxj,dyi,dyj,dzi,dzj;
	double kn,kb,kt;
	double g2lxx,g2lxy,g2lxz,g2lyx,g2lyy,g2lyz,g2lzx,g2lzy,g2lzz;
	double txi,txj,tyi,tyj,tzi,tzj;
	int p1, p2;

	for (size_t z = 0 ; z < nz_loc ; z++) {
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {

			int x_global = xsum + x;
			int y_global = ysum + y;
			int z_global = zsum + z;
				
		size_t k = ipx(x,y,z,nx_loc,ny_loc)/6;

		ip = ipx(x,y,z,nx_loc,ny_loc);
		ipy = ip + 1;
		ipz = ipy + 1;
		itx = ipz + 1;
		ity = itx + 1;
		itz = ity + 1;
								
		for (int l = 0 ; l < nodes[k].nb ; ++l) {
		
		jpx = nodes[k].ipx[l];

		jpy = jpx + 1;
		jpz = jpy + 1;
		jtx = jpz + 1;
		jty = jtx + 1;
		jtz = jty + 1;
			
		dxi = p[ip] - 	x_global * gridStep;
		dxj = p[jpx] - 	nodes[k].IX[l];
		dyi = p[ipy] - y_global * gridStep;
		dyj = p[jpy] - nodes[k].IY[l];
		dzi = p[ipz] - z_global * gridStep;
		dzj = p[jpz] - nodes[k].IZ[l];
			
		txi = p[itx];
		txj = p[jtx];
		tyi = p[ity];
		tyj = p[jty];
		tzi = p[itz];
		tzj = p[jtz];

		kn = _coef[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
		kb = _coefB[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
		kt = _coefT[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
			
		g2lxx = _G2Lxx[nodes[k].I[l]];
		g2lxy = _G2Lxy[nodes[k].I[l]];
		g2lxz = _G2Lxz[nodes[k].I[l]];
		g2lyx = _G2Lyx[nodes[k].I[l]];
		g2lyy = _G2Lyy[nodes[k].I[l]];
		g2lyz = _G2Lyz[nodes[k].I[l]];
		g2lzx = _G2Lzx[nodes[k].I[l]];
		g2lzy = _G2Lzy[nodes[k].I[l]];
		g2lzz = _G2Lzz[nodes[k].I[l]];
		
		df[ip] += 	0.5 * (2*dxi*(pow(g2lxx,2)*kn + (pow(g2lyx,2) + pow(g2lzx,2))*kt) - 
					2*dxj*(pow(g2lxx,2)*kn + (pow(g2lyx,2) + pow(g2lzx,2))*kt) + 
					2*dyi*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) - 
					2*dyj*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) + 
					2*dzi*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) - 
					2*dzj*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + 
					g2lzx*g2lzz*kt) + (-g2lyy*g2lzx + 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*tyi + (-g2lyy*g2lzx + 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*tyj + (-g2lyz*g2lzx + 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*tzi + (-g2lyz*g2lzx + g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*tzj);
		
		df[ipy] += 	0.5 * (2*dxi*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) - 
					2*dxj*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) + 
					2*dyi*(pow(g2lxy,2)*kn + (pow(g2lyy,2) + pow(g2lzy,2))*kt) - 
					2*dyj*(pow(g2lxy,2)*kn + (pow(g2lyy,2) + pow(g2lzy,2))*kt) + 
					2*dzi*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) - 
					2*dzj*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + 
					g2lzy*g2lzz*kt) + (g2lyy*g2lzx - 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*txi + (g2lyy*g2lzx - 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*txj + (-g2lyz*g2lzy + 
					g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tzi + (-g2lyz*g2lzy + g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tzj);
		
		df[ipz] += 	0.5 * (2*dxi*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) - 
					2*dxj*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) + 
					2*dyi*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) - 
					2*dyj*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) + 
					2*dzi*(pow(g2lxz,2)*kn + (pow(g2lyz,2) + pow(g2lzz,2))*kt) - 
					2*dzj*(pow(g2lxz,2)*kn + (pow(g2lyz,2) + pow(g2lzz,2))*kt) + (g2lyz*g2lzx - 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*txi + (g2lyz*g2lzx - 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*txj + (g2lyz*g2lzy - 
					g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tyi + (g2lyz*g2lzy - g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tyj);
		
		df[itx] += 	0.1667 * kt*_L0[nodes[k].I[l]]*(3*dyi*(g2lyy*g2lzx - g2lyx*g2lzy) + 
					3*dyj*(-g2lyy*g2lzx + g2lyx*g2lzy) + 
					3*dzi*(g2lyz*g2lzx - g2lyx*g2lzz) + 
					3*dzj*(-g2lyz*g2lzx + g2lyx*g2lzz) + 
					2*(pow(g2lyx,2) + pow(g2lzx,2))*_L0[nodes[k].I[l]]*txi + (pow(g2lyx,2) + pow(g2lzx,2))*_L0[nodes[k].I[l]]*txj + 
					2*(g2lyx*g2lyy + g2lzx*g2lzy)*_L0[nodes[k].I[l]]*tyi + (g2lyx*g2lyy + 
					g2lzx*g2lzy)*_L0[nodes[k].I[l]]*tyj + 
					2*(g2lyx*g2lyz + g2lzx*g2lzz)*_L0[nodes[k].I[l]]*tzi + (g2lyx*g2lyz + 
					g2lzx*g2lzz)*_L0[nodes[k].I[l]]*tzj);
		
		df[ity] += 	0.1667 * kt*_L0[nodes[k].I[l]]*(3*dxj*(g2lyy*g2lzx - g2lyx*g2lzy) + 
					3*dxi*(-g2lyy*g2lzx + g2lyx*g2lzy) + 
					3*dzi*(g2lyz*g2lzy - g2lyy*g2lzz) + 
					3*dzj*(-g2lyz*g2lzy + g2lyy*g2lzz) + 
					2*(g2lyx*g2lyy + g2lzx*g2lzy)*_L0[nodes[k].I[l]]*txi + (g2lyx*g2lyy + 
					g2lzx*g2lzy)*_L0[nodes[k].I[l]]*txj + 
					2*(pow(g2lyy,2) + pow(g2lzy,2))*_L0[nodes[k].I[l]]*tyi + (pow(g2lyy,2) + pow(g2lzy,2))*_L0[nodes[k].I[l]]*tyj + 
					2*(g2lyy*g2lyz + g2lzy*g2lzz)*_L0[nodes[k].I[l]]*tzi + (g2lyy*g2lyz + 
					g2lzy*g2lzz)*_L0[nodes[k].I[l]]*tzj);
					
		df[itz] += 	0.1667 * kt*_L0[nodes[k].I[l]]*(3*dxj*(g2lyz*g2lzx - g2lyx*g2lzz) + 
					3*dxi*(-g2lyz*g2lzx + g2lyx*g2lzz) + 
					3*dyj*(g2lyz*g2lzy - g2lyy*g2lzz) + 
					3*dyi*(-g2lyz*g2lzy + g2lyy*g2lzz) + 
					2*(g2lyx*g2lyz + g2lzx*g2lzz)*_L0[nodes[k].I[l]]*txi + (g2lyx*g2lyz + 
					g2lzx*g2lzz)*_L0[nodes[k].I[l]]*txj + 
					2*(g2lyy*g2lyz + g2lzy*g2lzz)*_L0[nodes[k].I[l]]*tyi + (g2lyy*g2lyz + 
					g2lzy*g2lzz)*_L0[nodes[k].I[l]]*tyj + 
					2*(pow(g2lyz,2) + pow(g2lzz,2))*_L0[nodes[k].I[l]]*tzi + (pow(g2lyz,2) + pow(g2lzz,2))*_L0[nodes[k].I[l]]*tzj);
						
		}
		
		for (int l = HALF ; l < HALF + nodes[k].nbo ; ++l) {
		
		jpx = nodes[k].ipx[l];		

		jpy = jpx + 1;
		jpz = jpy + 1;
		jtx = jpz + 1;
		jty = jtx + 1;
		jtz = jty + 1;
		
		dxi = p[ip] - 	x_global * gridStep;
		dxj = p[jpx] - 	nodes[k].IX[l];
		dyi = p[ipy] - y_global * gridStep;
		dyj = p[jpy] - nodes[k].IY[l];
		dzi = p[ipz] - z_global * gridStep;
		dzj = p[jpz] - nodes[k].IZ[l];
			
		txi = p[itx];
		txj = p[jtx];
		tyi = p[ity];
		tyj = p[jty];
		tzi = p[itz];
		tzj = p[jtz];

		kn = _coef[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
		kb = _coefB[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
		kt = _coefT[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
			
		g2lxx = _G2Lxx[nodes[k].I[l]];
		g2lxy = _G2Lxy[nodes[k].I[l]];
		g2lxz = _G2Lxz[nodes[k].I[l]];
		g2lyx = _G2Lyx[nodes[k].I[l]];
		g2lyy = _G2Lyy[nodes[k].I[l]];
		g2lyz = _G2Lyz[nodes[k].I[l]];
		g2lzx = _G2Lzx[nodes[k].I[l]];
		g2lzy = _G2Lzy[nodes[k].I[l]];
		g2lzz = _G2Lzz[nodes[k].I[l]];
		
		df[ip] += 	0.5 * (2*dxi*(pow(g2lxx,2)*kn + (pow(g2lyx,2) + pow(g2lzx,2))*kt) - 
					2*dxj*(pow(g2lxx,2)*kn + (pow(g2lyx,2) + pow(g2lzx,2))*kt) + 
					2*dyi*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) - 
					2*dyj*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) + 
					2*dzi*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) - 
					2*dzj*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + 
					g2lzx*g2lzz*kt) + (-g2lyy*g2lzx + 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*tyi + (-g2lyy*g2lzx + 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*tyj + (-g2lyz*g2lzx + 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*tzi + (-g2lyz*g2lzx + g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*tzj);
		
		df[ipy] += 	0.5 * (2*dxi*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) - 
					2*dxj*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) + 
					2*dyi*(pow(g2lxy,2)*kn + (pow(g2lyy,2) + pow(g2lzy,2))*kt) - 
					2*dyj*(pow(g2lxy,2)*kn + (pow(g2lyy,2) + pow(g2lzy,2))*kt) + 
					2*dzi*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) - 
					2*dzj*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + 
					g2lzy*g2lzz*kt) + (g2lyy*g2lzx - 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*txi + (g2lyy*g2lzx - 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*txj + (-g2lyz*g2lzy + 
					g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tzi + (-g2lyz*g2lzy + g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tzj);
		
		df[ipz] += 	0.5 * (2*dxi*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) - 
					2*dxj*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) + 
					2*dyi*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) - 
					2*dyj*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) + 
					2*dzi*(pow(g2lxz,2)*kn + (pow(g2lyz,2) + pow(g2lzz,2))*kt) - 
					2*dzj*(pow(g2lxz,2)*kn + (pow(g2lyz,2) + pow(g2lzz,2))*kt) + (g2lyz*g2lzx - 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*txi + (g2lyz*g2lzx - 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*txj + (g2lyz*g2lzy - 
					g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tyi + (g2lyz*g2lzy - g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tyj);
		
		df[itx] += 	0.1667 * kt*_L0[nodes[k].I[l]]*(3*dyi*(g2lyy*g2lzx - g2lyx*g2lzy) + 
					3*dyj*(-g2lyy*g2lzx + g2lyx*g2lzy) + 
					3*dzi*(g2lyz*g2lzx - g2lyx*g2lzz) + 
					3*dzj*(-g2lyz*g2lzx + g2lyx*g2lzz) + 
					2*(pow(g2lyx,2) + pow(g2lzx,2))*_L0[nodes[k].I[l]]*txi + (pow(g2lyx,2) + pow(g2lzx,2))*_L0[nodes[k].I[l]]*txj + 
					2*(g2lyx*g2lyy + g2lzx*g2lzy)*_L0[nodes[k].I[l]]*tyi + (g2lyx*g2lyy + 
					g2lzx*g2lzy)*_L0[nodes[k].I[l]]*tyj + 
					2*(g2lyx*g2lyz + g2lzx*g2lzz)*_L0[nodes[k].I[l]]*tzi + (g2lyx*g2lyz + 
					g2lzx*g2lzz)*_L0[nodes[k].I[l]]*tzj);
		
		df[ity] += 	0.1667 * kt*_L0[nodes[k].I[l]]*(3*dxj*(g2lyy*g2lzx - g2lyx*g2lzy) + 
					3*dxi*(-g2lyy*g2lzx + g2lyx*g2lzy) + 
					3*dzi*(g2lyz*g2lzy - g2lyy*g2lzz) + 
					3*dzj*(-g2lyz*g2lzy + g2lyy*g2lzz) + 
					2*(g2lyx*g2lyy + g2lzx*g2lzy)*_L0[nodes[k].I[l]]*txi + (g2lyx*g2lyy + 
					g2lzx*g2lzy)*_L0[nodes[k].I[l]]*txj + 
					2*(pow(g2lyy,2) + pow(g2lzy,2))*_L0[nodes[k].I[l]]*tyi + (pow(g2lyy,2) + pow(g2lzy,2))*_L0[nodes[k].I[l]]*tyj + 
					2*(g2lyy*g2lyz + g2lzy*g2lzz)*_L0[nodes[k].I[l]]*tzi + (g2lyy*g2lyz + 
					g2lzy*g2lzz)*_L0[nodes[k].I[l]]*tzj);
					
		df[itz] += 	0.1667 * kt*_L0[nodes[k].I[l]]*(3*dxj*(g2lyz*g2lzx - g2lyx*g2lzz) + 
					3*dxi*(-g2lyz*g2lzx + g2lyx*g2lzz) + 
					3*dyj*(g2lyz*g2lzy - g2lyy*g2lzz) + 
					3*dyi*(-g2lyz*g2lzy + g2lyy*g2lzz) + 
					2*(g2lyx*g2lyz + g2lzx*g2lzz)*_L0[nodes[k].I[l]]*txi + (g2lyx*g2lyz + 
					g2lzx*g2lzz)*_L0[nodes[k].I[l]]*txj + 
					2*(g2lyy*g2lyz + g2lzy*g2lzz)*_L0[nodes[k].I[l]]*tyi + (g2lyy*g2lyz + 
					g2lzy*g2lzz)*_L0[nodes[k].I[l]]*tyj + 
					2*(pow(g2lyz,2) + pow(g2lzz,2))*_L0[nodes[k].I[l]]*tzi + (pow(g2lyz,2) + pow(g2lzz,2))*_L0[nodes[k].I[l]]*tzj);

				}
			}
		}
	}
	
	for (size_t i = 0 ; i < ipImp.size() ; i++) df[ ipImp[i] ] = 0.0;
	for (size_t i = 0 ; i < ipTr.size() ; i++) 	df[ ipTr[i] ] += nTr[i];

	MPI_Barrier( comm3d );	
	MPI_Sendrecv (&(df[0]),1,EW_face,NeighBor[W],flag,&(df[nz_loc*nx_loc*ny_loc*ndof]), (nz_loc * ny_loc * ndof), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(df[ndof * (nx_loc - 1)]),1,EW_face,NeighBor[E],flag,&(df[nz_loc*nx_loc*ny_loc*ndof + nz_loc * ny_loc * ndof]), nz_loc * ny_loc * ndof, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(df[0]),1,bktype,NeighBor[F],flag,&(df[nz_loc * nx_loc * ny_loc * ndof + 2 * nz_loc * ny_loc * ndof]), ((nx_loc + 2) * nz_loc * ndof), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(df[0]),1,frtype,NeighBor[B],flag,&(df[nz_loc * nx_loc * ny_loc * ndof + 2 * nz_loc * ny_loc * ndof + (nx_loc + 2) * nz_loc * ndof]), ((nx_loc + 2) * nz_loc * ndof), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(df[0]),1,tptype,NeighBor[N],flag,&(df[nx_loc * ny_loc * nz_loc * ndof + 2 * nz_loc * ny_loc * ndof + (nx_loc + 2) * nz_loc * ndof * 2]),(nx_fill * ny_fill * ndof),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(df[0]),1,bttype,NeighBor[S],flag,&(df[nx_loc * ny_loc * nz_loc * ndof + 2 * ny_loc * nz_loc * ndof + 2 * (nx_loc + 2) * nz_loc * ndof + nx_fill * ny_fill * ndof]), (nx_fill * ny_fill * ndof),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);
}


// Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and
// resets p to where the function func(p) takes on a minimum along the direction xi from p,
// and replaces xi by the actual vector displacement that p was moved. Also returns as fret
// the value of func at the returned location p. This is actually all accomplished by calling the
// routines mnbrak and brent.
void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
{
	int j;
	double xx, xmin, fx, fb, fa, bx, ax;
	ncom = n;                                     
	pcom = new double[n];
	xicom = new double[n];
	nrfunc = func;
	
	for (j = 0; j < n; j++)
	{
		pcom[j] = p[j];
		xicom[j] = xi[j];
	}
	
	ax = 0.0;                                    
	xx = 1.0;
	
	mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
	*fret = brent(ax, xx, bx, f1dim, TOL, &xmin);

	for (j = 0; j < n; j++)
	{
		xi[j] *= xmin;
		p[j] += xi[j];
	}

	delete[] xicom;
	xicom = 0;
	delete[] pcom;
	pcom = 0;
}

double f1dim(double x)
{
	int j;
	double f,f_tot,*xt;
	xt = new double[ncom];
	for (j = 0; j < ncom; j++) xt[j] = pcom[j] + x * xicom[j];
	
	f = (*nrfunc)(xt);
	MPI_Allreduce(&f, &f_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	
	delete[] xt;
	xt = 0;
	return f_tot;
}

// Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a
// function func, using its gradient as calculated by a routine dfunc. The convergence tolerance
// on the function value is input as ftol. Returned quantities are p (the location of the minimum),
// iter (the number of iterations that were performed), and fret (the minimum value of the
// function). The routine linmin is called to perform line minimizations.
// (p. 447)
void frprmn(double p[], size_t n, double ftol, size_t fmaxIter, int *iter, double *fret, double (*func)(double []), void (*dfunc)(double [], double []), int rrkk){
	const char* endcond = "iterations";
	size_t j, its;
	double gg, gam, fp,fpi, dgg;
	double *g, *h, *xi;
	
	g = new double[n];
	h = new double[n];
	xi = new double[n];
	
	fpi = (*func)(p);  
	(*dfunc)(p, xi);
		
	MPI_Allreduce(&fp, &fpi, 1, MPI_DOUBLE, MPI_SUM,comm3d);
		
	for (j = 0; j < n; j++){
		g[j] = -xi[j];
		xi[j] = h[j] = g[j];
	}
	
	for (its = 1; its <= fmaxIter; its++) { 
		*iter = its;
		linmin(p, xi, n, fret, func); 

		//if((its% 50)==0){
		if(rrkk == 0) cout << "iteration " << its << ", fret " << *fret << endl;
		//}
		
		if (2.0 * fabs(*fret - fp) <= ftol * (fabs(*fret) + fabs(fp) + EPS)){
			endcond = "tolerance";
			delete[] g;
			g = 0;
			delete[] h;
			h = 0;
			delete[] xi;
			xi = 0;
			break;
		}

	fp = *fret;
	(*dfunc)(p, xi);
				
		dgg = gg = 0.0;
		
		for (j = 0; j < (nx_loc*ny_loc*nz_loc*ndof); j++){
			gg += g[j] * g[j];
			//dgg += xi[j]*xi[j];          // This statement for Fletcher-Reeves.
			dgg += (xi[j] + g[j]) * xi[j];  // This statement for Polak-Ribiere.
		}
		
		double ggg = 0.0;
		double gdgg = 0.0;
		MPI_Barrier(comm3d);
		MPI_Allreduce(&gg, &ggg, 1, MPI_DOUBLE, MPI_SUM, comm3d);
		MPI_Allreduce(&dgg, &gdgg, 1, MPI_DOUBLE, MPI_SUM, comm3d);
				
		if (ggg == 0.0){
			endcond = "threshold";
			delete[] g;
			g = 0;
			delete[] h;
			h = 0;
			delete[] xi;
			xi = 0;    
			break;
		}
		
		gam = gdgg/ggg;
		
		for (j = 0; j < n; j++){
			g[j] = -xi[j];
			xi[j] = h[j] = g[j] + gam * h[j];
		}

	}
	
	return;
}

// Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
// between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
// the minimum to a fractional precision of about tol using Brent method. The abscissa of
// the minimum is returned as xmin, and the minimum function value is returned as brent, the
// returned function value.
double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin)
{
	int iter;
	double a, b, d = 0.0f, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm, ftmp, ftmp_tot;
	double e = 0.0;                                  

	a = (ax < cx ? ax : cx);                        	
	b = (ax > cx ? ax : cx);                        	
	x = w = v = bx;                                    	
	fw = fv = fx = (*f)(x);

	for (iter = 1; iter <= ITMAX_BRENT; iter++) {      	 
		xm = 0.5 * (a + b);
		tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
		if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) { 
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1) {                 
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0) p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			else {
				d = p / q;                       
				u = x + d;
				if (u - a < tol2 || b - u < tol2)
					d = SIGN(tol1, xm - x);
			}
		}
		else {
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		}
		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu = (*f)(u);
		
		if (fu <= fx) {                  
			if (u >= x) a = x;
			else b = x;
			SHFT(v, w, x, u)                   
			SHFT(fv, fw, fx, fu)
		}
		else {
			if (u < x) a = u;
			else b = u;
			if (fu <= fw || w == x) {
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}                                   
	}                                      
	cout << "Too many iterations in brent" << endl;
	*xmin = x;                                   
	return fx;
}

// Given a function func, and given distinct initial points ax and bx, this routine searches in
// the downhill direction (defined by the function as evaluated at the initial points) and returns
// new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
// values at the three points, fa, fb, and fc.
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double))
{
	double ulim, u, r, q, fu, dum;
	*fa = (*func)(*ax);
	*fb = (*func)(*bx);

	if (*fb > *fa) { 
		SHFT(dum, *ax, *bx, dum) 
		SHFT(dum, *fb, *fa, dum)
	}
	*cx = (*bx) + GOLD * (*bx - *ax);
	*fc = (*func)(*cx);

	while (*fb > *fc) { 
		r = (*bx - *ax) * (*fb - *fc);
		q = (*bx - *cx) * (*fb - *fa);
		u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) /
		    (2.0 * SIGN(FMAX(fabs(q - r), TINY), q - r));
		ulim = (*bx) + GLIMIT * (*cx - *bx);
		if ((*bx - u) * (u - *cx) > 0.0) { 
		fu = (*func)(u);
			
			if (fu < *fc) { 
				*ax = (*bx);
				*bx = u;
				*fa = (*fb);
				*fb = fu;
				return;
			}
			else if (fu > *fb) {
				*cx = u;
				*fc = fu;
				return;
			}
			u = (*cx) + GOLD * (*cx - *bx);  
			fu = (*func)(u);

		}
		else if ((*cx - u) * (u - ulim) > 0.0) {	
		fu = (*func)(u);
			
			if (fu < *fc) {
				SHFT(*bx, *cx, u, *cx + GOLD * (*cx - *bx))
				SHFT(*fb, *fc, fu, (*func)(u))
			}
		}
		else if ((u - ulim) * (ulim - *cx) >= 0.0) {	
			u = ulim; 
			fu = (*func)(u);
		}
		else {	
			u = (*cx) + GOLD * (*cx - *bx);
			fu = (*func)(u);
		}
		SHFT(*ax, *bx, *cx, u) 
		SHFT(*fa, *fb, *fc, fu)
	}
}


void cylinder(int idincl, double ideta, double xc, double yc, double zc, double rad, size_t nx_loc, size_t ny_loc, size_t nz_loc, int xsum, int ysum, int zsum, int xshift, int yshift, int zshift){

    double crit = 0;

	xc = xc + xshift;
	yc = yc + yshift;
	zc = zc + zshift;
    
    	if(rad > 0){
 	for (size_t z = 0 ; z < nz_loc ; z++){					
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);

				x_global = x_global - 1;
				y_global = y_global;
				z_global = z_global;

				crit=(1.*z_global*gridStep-1.*zc)*(1.*z_global*gridStep-1.*zc)+(1.*y_global*gridStep-1.*yc)*(1.*y_global*gridStep-1.*yc);
				
				if(crit<=1.*rad*rad){
   				Id[ipx(x,y,z,nx_loc,ny_loc)/6] = idincl;
				Eta[ipx(x,y,z,nx_loc,ny_loc)/6] = ideta;
    				}
   			 }
    		}
   	 } 
    }
		

	MPI_Barrier(comm3d);
	MPI_Sendrecv (&(Id[0]),1,EW_faceID,NeighBor[W],flag,&(Id[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_INT,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[1 * (nx_loc - 1)]),1,EW_faceID,NeighBor[E],flag,&(Id[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_INT,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(Id[0]),1,bktypeID,NeighBor[F],flag,&(Id[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_INT,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,frtypeID,NeighBor[B],flag,&(Id[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_INT,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,tptypeID,NeighBor[N],flag,&(Id[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_INT,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,bttypeID,NeighBor[S],flag,&(Id[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_INT,NeighBor[N],flag,comm3d,&status);

	MPI_Sendrecv (&(Eta[0]),1,EW_faceETA,NeighBor[W],flag,&(Eta[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[1 * (nx_loc - 1)]),1,EW_faceETA,NeighBor[E],flag,&(Eta[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(Eta[0]),1,bktypeETA,NeighBor[F],flag,&(Eta[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,frtypeETA,NeighBor[B],flag,&(Eta[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,tptypeETA,NeighBor[N],flag,&(Eta[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,bttypeETA,NeighBor[S],flag,&(Eta[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);
	
}

void  Import_Structure(const string& _name, int idincl, double ideta, size_t npart,size_t nx_loc, size_t ny_loc, size_t nz_loc, int xsum, int ysum, int zsum, int xshift, int yshift, int zshift){
	
	double xc,yc,zc,crit,rad;

	fstream file(_name);

	for (size_t k = 0 ; k < npart ; k++){
	file >> xc >> yc >> zc >> rad;
		
	xc = xc + xshift;
	yc = yc + yshift;
	zc = zc + zshift;
		
 	for (size_t z = 0 ; z < nz_loc ; z++){					
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = xsum + x;
				int y_global = ysum + y;
				int z_global = zsum + z;

				x_global = x_global - 1;
				y_global = y_global - 1;
				z_global = z_global - 1;

				crit=(1.*x_global*gridStep-1.*xc)*(1.*x_global*gridStep-1.*xc)+(1.*y_global*gridStep-1.*yc)*(1.*y_global*gridStep-1.*yc)+(1.*z_global*gridStep-1.*zc)*(1.*z_global*gridStep-1.*zc);

				if(crit<=1.*rad*rad){
   				Id[ipx(x,y,z,nx_loc,ny_loc)/6] = idincl;
				Eta[ipx(x,y,z,nx_loc,ny_loc)/6] = ideta;
    				}
    			}
    		}
    	}
	}
			
	MPI_Barrier(comm3d);
	MPI_Sendrecv (&(Id[0]),1,EW_faceID,NeighBor[W],flag,&(Id[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_INT,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[1 * (nx_loc - 1)]),1,EW_faceID,NeighBor[E],flag,&(Id[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_INT,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(Id[0]),1,bktypeID,NeighBor[F],flag,&(Id[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_INT,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,frtypeID,NeighBor[B],flag,&(Id[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_INT,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,tptypeID,NeighBor[N],flag,&(Id[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_INT,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,bttypeID,NeighBor[S],flag,&(Id[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_INT,NeighBor[N],flag,comm3d,&status);

	MPI_Sendrecv (&(Eta[0]),1,EW_faceETA,NeighBor[W],flag,&(Eta[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[1 * (nx_loc - 1)]),1,EW_faceETA,NeighBor[E],flag,&(Eta[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(Eta[0]),1,bktypeETA,NeighBor[F],flag,&(Eta[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,frtypeETA,NeighBor[B],flag,&(Eta[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,tptypeETA,NeighBor[N],flag,&(Eta[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,bttypeETA,NeighBor[S],flag,&(Eta[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);	
}

/*
void disorder(int idincl, double ideta, double xc, double yc, double zc, double rad, size_t nx_loc, size_t ny_loc, size_t nz_loc, int xsum, int ysum, int zsum){

    double crit = 0;
    
    if(rad > 0){
 	for (size_t z = 0 ; z < nz_loc ; z++){					
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);

				crit=(1.*x_global*gridStep-1.*xc)*(1.*x_global*gridStep-1.*xc)+(1.*y_global*gridStep-1.*yc)*(1.*y_global*gridStep-1.*yc)+(1.*z_global*gridStep-1.*zc)*(1.*z_global*gridStep-1.*zc);

				if(crit<=1.*rad*rad){
   				Id[ipx(x,y,z,nx_loc,ny_loc)/6] = idincl;
				Eta[ipx(x,y,z,nx_loc,ny_loc)/6] = ideta;
    				}
    			}
    		}
    	} 
    }
		

	MPI_Barrier(comm3d);
	MPI_Sendrecv (&(Id[0]),1,EW_faceID,NeighBor[W],flag,&(Id[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_INT,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[1 * (nx_loc - 1)]),1,EW_faceID,NeighBor[E],flag,&(Id[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_INT,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(Id[0]),1,bktypeID,NeighBor[F],flag,&(Id[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_INT,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,frtypeID,NeighBor[B],flag,&(Id[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_INT,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,tptypeID,NeighBor[N],flag,&(Id[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_INT,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,bttypeID,NeighBor[S],flag,&(Id[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_INT,NeighBor[N],flag,comm3d,&status);

	MPI_Sendrecv (&(Eta[0]),1,EW_faceETA,NeighBor[W],flag,&(Eta[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[1 * (nx_loc - 1)]),1,EW_faceETA,NeighBor[E],flag,&(Eta[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(Eta[0]),1,bktypeETA,NeighBor[F],flag,&(Eta[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,frtypeETA,NeighBor[B],flag,&(Eta[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,tptypeETA,NeighBor[N],flag,&(Eta[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,bttypeETA,NeighBor[S],flag,&(Eta[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);
	
}
*/


/*
void  disorder(int idincl, double ideta, double xc, double yc, double zc, double rad, size_t nx_loc, size_t ny_loc, size_t nz_loc, int xsum, int ysum, int zsum, int xshift, int yshift, int zshift){
	
	double crit;
		
	xc = xc + xshift;
	yc = yc + yshift;
	zc = zc + zshift;
		
 	for (size_t z = 0 ; z < nz_loc ; z++){					
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = xsum + x;
				int y_global = ysum + y;
				int z_global = zsum + z;

				x_global = x_global - 1;
				y_global = y_global - 1;
				z_global = z_global - 1;

				crit=(1.*x_global*gridStep-1.*xc)*(1.*x_global*gridStep-1.*xc)+(1.*y_global*gridStep-1.*yc)*(1.*y_global*gridStep-1.*yc)+(1.*z_global*gridStep-1.*zc)*(1.*z_global*gridStep-1.*zc);

				if(crit<=1.*rad*rad){
   				Id[ipx(x,y,z,nx_loc,ny_loc)/6] = idincl;
				Eta[ipx(x,y,z,nx_loc,ny_loc)/6] = ideta;
    				}
    			}
    		}
    	}
				
	MPI_Barrier(comm3d);
	MPI_Sendrecv (&(Id[0]),1,EW_faceID,NeighBor[W],flag,&(Id[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_INT,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[1 * (nx_loc - 1)]),1,EW_faceID,NeighBor[E],flag,&(Id[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_INT,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(Id[0]),1,bktypeID,NeighBor[F],flag,&(Id[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_INT,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,frtypeID,NeighBor[B],flag,&(Id[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_INT,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,tptypeID,NeighBor[N],flag,&(Id[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_INT,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,bttypeID,NeighBor[S],flag,&(Id[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_INT,NeighBor[N],flag,comm3d,&status);

	MPI_Sendrecv (&(Eta[0]),1,EW_faceETA,NeighBor[W],flag,&(Eta[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[1 * (nx_loc - 1)]),1,EW_faceETA,NeighBor[E],flag,&(Eta[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(Eta[0]),1,bktypeETA,NeighBor[F],flag,&(Eta[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,frtypeETA,NeighBor[B],flag,&(Eta[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,tptypeETA,NeighBor[N],flag,&(Eta[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,bttypeETA,NeighBor[S],flag,&(Eta[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);	
}
*/


void computeStress(size_t k,size_t x, size_t y, size_t z, tensor & stress)
{
	size_t ip = k * 6;
	size_t ipy = ip + 1;
	size_t ipz = ipy + 1;
	size_t itx = ipz + 1;
	size_t ity = itx + 1;
	size_t itz = ity + 1;
	
	size_t jpx, jpy, jpz, jtx, jty, jtz;
	double txi,txj,tyi,tyj,tzi,tzj;
	double dxi,dxj,dyi,dyj,dzi,dzj;
	double kn,kb,kt;
	double g2lxx,g2lxy,g2lxz,g2lyx,g2lyy,g2lyz,g2lzx,g2lzy,g2lzz;
	double fx, fy, fz;
	
	int x_global = xsum + x;
	int y_global = ysum + y;
	int z_global = zsum + z;

	stress.xx=stress.xy=stress.xz=stress.yx=stress.yy=stress.yz=stress.zx=stress.zy=stress.zz=0.0;
	
	for (int l = 0 ; l < nodes[k].nb ; ++l) {
			
	jpx = nodes[k].ipx[l];

	jpy = jpx + 1;
	jpz = jpy + 1;
	jtx = jpz + 1;
	jty = jtx + 1;
	jtz = jty + 1;
			
	dxi = p[ip] -  x_global * gridStep;
	dxj = p[jpx] - nodes[k].IX[l];
	dyi = p[ipy] - y_global * gridStep;
	dyj = p[jpy] - nodes[k].IY[l];
	dzi = p[ipz] - z_global * gridStep;
	dzj = p[jpz] - nodes[k].IZ[l];				

	txi = p[itx];
	txj = p[jtx];
	tyi = p[ity];
	tyj = p[jty];
	tzi = p[itz];
	tzj = p[jtz];

	kn = _coef[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
	kb = _coefB[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
	kt = _coefT[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
			
	g2lxx = _G2Lxx[nodes[k].I[l]];
	g2lxy = _G2Lxy[nodes[k].I[l]];
	g2lxz = _G2Lxz[nodes[k].I[l]];
	g2lyx = _G2Lyx[nodes[k].I[l]];
	g2lyy = _G2Lyy[nodes[k].I[l]];
	g2lyz = _G2Lyz[nodes[k].I[l]];
	g2lzx = _G2Lzx[nodes[k].I[l]];
	g2lzy = _G2Lzy[nodes[k].I[l]];
	g2lzz = _G2Lzz[nodes[k].I[l]];
		
		fx = 	0.5 * (2*dxi*(pow(g2lxx,2)*kn + (pow(g2lyx,2) + pow(g2lzx,2))*kt) - 
					2*dxj*(pow(g2lxx,2)*kn + (pow(g2lyx,2) + pow(g2lzx,2))*kt) + 
					2*dyi*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) - 
					2*dyj*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) + 
					2*dzi*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) - 
					2*dzj*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + 
					g2lzx*g2lzz*kt) + (-g2lyy*g2lzx + 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*tyi + (-g2lyy*g2lzx + 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*tyj + (-g2lyz*g2lzx + 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*tzi + (-g2lyz*g2lzx + g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*tzj)
					+ nodes[ipx(x,y,z,nx_loc,ny_loc)/6].FX[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]];
		
		fy = 	0.5 * (2*dxi*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) - 
					2*dxj*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) + 
					2*dyi*(pow(g2lxy,2)*kn + (pow(g2lyy,2) + pow(g2lzy,2))*kt) - 
					2*dyj*(pow(g2lxy,2)*kn + (pow(g2lyy,2) + pow(g2lzy,2))*kt) + 
					2*dzi*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) - 
					2*dzj*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + 
					g2lzy*g2lzz*kt) + (g2lyy*g2lzx - 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*txi + (g2lyy*g2lzx - 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*txj + (-g2lyz*g2lzy + 
					g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tzi + (-g2lyz*g2lzy + g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tzj)
					+ nodes[ipx(x,y,z,nx_loc,ny_loc)/6].FY[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]];
		
		fz = 	0.5 * (2*dxi*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) - 
					2*dxj*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) + 
					2*dyi*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) - 
					2*dyj*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) + 
					2*dzi*(pow(g2lxz,2)*kn + (pow(g2lyz,2) + pow(g2lzz,2))*kt) - 
					2*dzj*(pow(g2lxz,2)*kn + (pow(g2lyz,2) + pow(g2lzz,2))*kt) + (g2lyz*g2lzx - 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*txi + (g2lyz*g2lzx - 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*txj + (g2lyz*g2lzy - 
					g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tyi + (g2lyz*g2lzy - g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tyj)
					+ nodes[ipx(x,y,z,nx_loc,ny_loc)/6].FZ[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]];
										        		
    stress.xx -= fx * pow(gridStep,2) * _nx[nodes[k].I[l]];
    stress.xy -= fx * pow(gridStep,2) * _ny[nodes[k].I[l]];
    stress.xz -= fx * pow(gridStep,2) * _nz[nodes[k].I[l]];
    stress.yx -= fy * pow(gridStep,2) * _nx[nodes[k].I[l]];
    stress.yy -= fy * pow(gridStep,2) * _ny[nodes[k].I[l]];
    stress.yz -= fy * pow(gridStep,2) * _nz[nodes[k].I[l]];
    stress.zx -= fz * pow(gridStep,2) * _nx[nodes[k].I[l]];
    stress.zy -= fz * pow(gridStep,2) * _ny[nodes[k].I[l]];
    stress.zz -= fz * pow(gridStep,2) * _nz[nodes[k].I[l]];

	}
	
	for (int l = HALF ; l < HALF + nodes[k].nbo ; ++l) {
		
	jpx = nodes[k].ipx[l];

	jpy = jpx + 1;
	jpz = jpy + 1;
	jtx = jpz + 1;
	jty = jtx + 1;
	jtz = jty + 1;
			
	dxi = p[ip] -  x_global * gridStep;
	dxj = p[jpx] - nodes[k].IX[l];
	dyi = p[ipy] - y_global * gridStep;
	dyj = p[jpy] - nodes[k].IY[l];
	dzi = p[ipz] - z_global * gridStep;
	dzj = p[jpz] - nodes[k].IZ[l];				

	txi = p[itx];
	txj = p[jtx];
	tyi = p[ity];
	tyj = p[jty];
	tzi = p[itz];
	tzj = p[jtz];
  
	kn = _coef[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
	kb = _coefB[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
	kt = _coefT[nodes[k].I[l]] * Ktable[Id[k]][Id[jpx/6]];
			
	g2lxx = _G2Lxx[nodes[k].I[l]];
	g2lxy = _G2Lxy[nodes[k].I[l]];
	g2lxz = _G2Lxz[nodes[k].I[l]];
	g2lyx = _G2Lyx[nodes[k].I[l]];
	g2lyy = _G2Lyy[nodes[k].I[l]];
	g2lyz = _G2Lyz[nodes[k].I[l]];
	g2lzx = _G2Lzx[nodes[k].I[l]];
	g2lzy = _G2Lzy[nodes[k].I[l]];
	g2lzz = _G2Lzz[nodes[k].I[l]];
		
		fx = 	0.5 * (2*dxi*(pow(g2lxx,2)*kn + (pow(g2lyx,2) + pow(g2lzx,2))*kt) - 
					2*dxj*(pow(g2lxx,2)*kn + (pow(g2lyx,2) + pow(g2lzx,2))*kt) + 
					2*dyi*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) - 
					2*dyj*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) + 
					2*dzi*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) - 
					2*dzj*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + 
					g2lzx*g2lzz*kt) + (-g2lyy*g2lzx + 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*tyi + (-g2lyy*g2lzx + 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*tyj + (-g2lyz*g2lzx + 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*tzi + (-g2lyz*g2lzx + g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*tzj)
					+ nodes[ipx(x,y,z,nx_loc,ny_loc)/6].FX[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]];
		
		fy = 	0.5 * (2*dxi*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) - 
					2*dxj*(g2lxx*g2lxy*kn + g2lyx*g2lyy*kt + g2lzx*g2lzy*kt) + 
					2*dyi*(pow(g2lxy,2)*kn + (pow(g2lyy,2) + pow(g2lzy,2))*kt) - 
					2*dyj*(pow(g2lxy,2)*kn + (pow(g2lyy,2) + pow(g2lzy,2))*kt) + 
					2*dzi*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) - 
					2*dzj*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + 
					g2lzy*g2lzz*kt) + (g2lyy*g2lzx - 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*txi + (g2lyy*g2lzx - 
					g2lyx*g2lzy)*kt*_L0[nodes[k].I[l]]*txj + (-g2lyz*g2lzy + 
					g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tzi + (-g2lyz*g2lzy + g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tzj)
					+ nodes[ipx(x,y,z,nx_loc,ny_loc)/6].FY[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]];
		
		fz = 	0.5 * (2*dxi*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) - 
					2*dxj*(g2lxx*g2lxz*kn + g2lyx*g2lyz*kt + g2lzx*g2lzz*kt) + 
					2*dyi*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) - 
					2*dyj*(g2lxy*g2lxz*kn + g2lyy*g2lyz*kt + g2lzy*g2lzz*kt) + 
					2*dzi*(pow(g2lxz,2)*kn + (pow(g2lyz,2) + pow(g2lzz,2))*kt) - 
					2*dzj*(pow(g2lxz,2)*kn + (pow(g2lyz,2) + pow(g2lzz,2))*kt) + (g2lyz*g2lzx - 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*txi + (g2lyz*g2lzx - 
					g2lyx*g2lzz)*kt*_L0[nodes[k].I[l]]*txj + (g2lyz*g2lzy - 
					g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tyi + (g2lyz*g2lzy - g2lyy*g2lzz)*kt*_L0[nodes[k].I[l]]*tyj)
					+ nodes[ipx(x,y,z,nx_loc,ny_loc)/6].FZ[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]];
										        		
    stress.xx -= fx * pow(gridStep,2) * _nx[nodes[k].I[l]];
    stress.xy -= fx * pow(gridStep,2) * _ny[nodes[k].I[l]];
    stress.xz -= fx * pow(gridStep,2) * _nz[nodes[k].I[l]];
    stress.yx -= fy * pow(gridStep,2) * _nx[nodes[k].I[l]];
    stress.yy -= fy * pow(gridStep,2) * _ny[nodes[k].I[l]];
    stress.yz -= fy * pow(gridStep,2) * _nz[nodes[k].I[l]];
    stress.zx -= fz * pow(gridStep,2) * _nx[nodes[k].I[l]];
    stress.zy -= fz * pow(gridStep,2) * _ny[nodes[k].I[l]];
    stress.zz -= fz * pow(gridStep,2) * _nz[nodes[k].I[l]];

	}
		
	double V = 2.* gridStep * gridStep * gridStep;
	stress.xx /= V;
	stress.xy /= V;
	stress.xz /= V;
	stress.yx /= V;
	stress.yy /= V;
	stress.yz /= V;
	stress.zx /= V;
	stress.zy /= V;
	stress.zz /= V;
}

void force_init(){
	
	 for (size_t z = 0 ; z < nz_loc ; z++){					
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {
				for (int l = 0 ; l < Q ; l++) {
					nodes[ipx(x,y,z,nx_loc,ny_loc)/6].FX[l] = 0.0;
					nodes[ipx(x,y,z,nx_loc,ny_loc)/6].FY[l] = 0.0;
					nodes[ipx(x,y,z,nx_loc,ny_loc)/6].FZ[l] = 0.0;
				}	
			}
		}
	 }
}

void FabricTensor(int idp, int rrkk)
{
	
	int ph;
	double Fa11, Fa22, Fa33;
	Fa11 = Fa22 = Fa33 = 0.0;

	for (size_t z = 0 ; z < nz_loc ; z++) {					
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	ph = Id[ipx(x,y,z,nx_loc,ny_loc)/6];

	if(ph == idp){
	for (int l = 0 ; l < nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nb ; ++l) {
	Fa11+= (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
	Fa22+= (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
	Fa33+= (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
	}
	for (int l = HALF ; l < HALF + nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nbo ; ++l) {	
	Fa11+= (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
	Fa22+= (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
	Fa33+= (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
	}
	}

	}
	}
	}

MPI_Barrier(comm3d);
MPI_Allreduce(&Fa11, &f11, 1, MPI_DOUBLE, MPI_SUM, comm3d);
MPI_Allreduce(&Fa22, &f22, 1, MPI_DOUBLE, MPI_SUM, comm3d);
MPI_Allreduce(&Fa33, &f33, 1, MPI_DOUBLE, MPI_SUM, comm3d);

}

void ApplyPorePressure(int idp, int ids, double Pressure, double numP, int rrkk){
	
	double Force = -Pressure * 6 * (numP) / (f11 + f22 + f33);
	int ph;

 	for (size_t z = 0 ; z < nz_loc ; z++){					
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {
				
		int x_global = xsum + x;
		int y_global = ysum + y;
		int z_global = zsum + z;
				
		ph = Id[ipx(x,y,z,nx_loc,ny_loc)/6];
						
		if (ph == ids){
			
		for (int l = 0 ; l < nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nb ; ++l) {

		if (Id[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l]/6] == idp){
		force[ipx(x,y,z,nx_loc,ny_loc)] += (Force/_L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
		force[ipx(x,y,z,nx_loc,ny_loc) + 1] += (Force/_L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
		force[ipx(x,y,z,nx_loc,ny_loc) + 2] += (Force/_L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
		}
		}
		for (int l = HALF ; l < HALF + nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nbo ; ++l) {

		if (Id[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l]/6] == idp){
		force[ipx(x,y,z,nx_loc,ny_loc)] += (Force/_L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
		force[ipx(x,y,z,nx_loc,ny_loc) + 1] += (Force/_L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
		force[ipx(x,y,z,nx_loc,ny_loc) + 2] += (Force/_L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) * (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);		
		}
		}
		}

		}
		}
	}
}


void Force_BC(int id0, int id1, int rrkk,size_t nx_loc, size_t ny_loc, size_t nz_loc, int xsum, int ysum, int zsum){

	vector<size_t> IntIndList;
		
 	for (size_t z = 0 ; z < nz_loc ; z++){					
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {
				
		int x_global = xsum + x;
		int y_global = ysum + y;
		int z_global = zsum + z;

		if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == id1){
		for (int l = 0 ; l < nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nb ; ++l) {
		if	(Id[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l]/6] == id0){
		IntIndList.push_back(ipx(x,y,z,nx_loc,ny_loc));
		}
		}
		for (int l = HALF ; l < HALF + nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nbo ; ++l) {
		if	(Id[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l]/6] == id0){
		IntIndList.push_back(ipx(x,y,z,nx_loc,ny_loc));
		}
		}
		}
			
		}
		}
		}
		
		
sort( IntIndList.begin(), IntIndList.end() );
IntIndList.erase( unique( IntIndList.begin(), IntIndList.end() ), IntIndList.end() );
	
for (size_t i = 0 ; i < IntIndList.size() ; i++){

nTr.push_back(force[IntIndList[i]]);
nTr.push_back(force[IntIndList[i] + 1]);
nTr.push_back(force[IntIndList[i] + 2]);
	
ipTr.push_back(IntIndList[i]);
ipTr.push_back(IntIndList[i] + 1);
ipTr.push_back(IntIndList[i] + 2);
	
}
}

double vol_frac(double tmp[], int idp, int rrkk, int xls, int xrs, int yls, int yrs, int zds, int zus){

	double count = 0.0;
	double count_A = 0.0;
	double count_B = 0.0;
	bool res_cond, box_cond;

 	for (size_t z = 0 ; z < nz_loc ; z++) {					
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);

          res_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) && ( z_global >= zds ) && ( z_global <  (nz) - zus) ? true: false;//if true -> continue 
          box_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) && ( z_global >= zds ) && ( z_global <  (nz) - zus) ? false: true;//if false -> continue 

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == idp){
	count++;
	}

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == idp && res_cond){
	count_A++;
	}

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == idp && box_cond){
	count_B++;
	}

	}
	}
	}

	double fs,fp, count_tot, count_tot_A, count_tot_B;
	MPI_Barrier(comm3d);
	MPI_Allreduce(&count, &count_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&count_A, &count_tot_A, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&count_B, &count_tot_B, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	double tot = nx * ny * nz;
	double den = pow((cbrt(nx * ny * nz)-1),3);
	fp = (count_tot)/den;
	fs = 1.0 - fp;

	tmp[0] = count_tot;
	tmp[1] = count_tot_A;
	tmp[2] = count_tot_B;
	return fp;
}


void Elastic_Input_ISO(double M, double NU, double N){
	
	double YOUNG = (1 - pow(NU,2)) * M;
	double SHEAR = YOUNG / 2.0 / (1.0 + NU);
	double BULK = YOUNG / 3.0 / (1 - 2.0*NU);
	double LAME = (3.0 * BULK * NU) / (1.0 + NU);
	double C11 = LAME + 2.0 * SHEAR;
	double C12 = LAME;
	double C13 = LAME;
	double C33 = C11;
	double C44 = SHEAR;

	double k1 = (1.0/M/pow((1 + N),2)) * pow(N,2) * (C33 - 2.0 * C12);
	double k2 = (1.0/M/(1 + N)) * N * (C12);
	double k4 = (1.0/M/pow((1+N),2)) * pow(N,2) * (-1.5 * C12 + 0.5 * C33);

	for(int i = 0 ; i < HALF ; i++){
	if(i < 3){
	_coef[i] = k1;
	_coefB[i] = k4;
	_coefT[i] = k4;
	
	_coef[i + HALF] = k1;
	_coefB[i + HALF] = k4;
	_coefT[i + HALF] = k4;
	}
	if(i > 2 && i < 9){
	_coef[i] = k2;
	_coef[i + HALF] = k2;
	}
	if(i > 8 && i < 13){
	_coef[i] = 0.0;
	_coef[i + HALF] = 0.0;
	}
	if(i > 2 && i < 13){
	_coefB[i] = 0.0;
	_coefT[i] = 0.0;
	_coefB[i + HALF] = 0.0;
	_coefT[i + HALF] = 0.0;
	}
	}
	
}


void save_MPI_info(int rrkk,int nx_loc,int ny_loc, int nz_loc, int xcoor,int ycoor,int zcoor,int xs,int xsp,int ys,int ysp,int zs,int zsp){

    FILE * sortie;
    char nomfic[256];
	
    sprintf(nomfic, "mpi_info.%04i", rrkk);
    sortie = fopen(nomfic, "w+");
    fprintf(sortie,"%i %i %i %i %i %i %i %i %i %i %i %i %i\n",rrkk,nx_loc,ny_loc,nz_loc,xcoor,ycoor,zcoor,xs,xs+xsp,ys,ys+ysp,zs,zs+zsp);
    fclose(sortie);	

}

void save_output(int rrkk){
	
    FILE * sortie;
    char nomfic[256];
	
    sprintf(nomfic, "out.%04i", rrkk);
    sortie = fopen(nomfic, "w+");

    double u_x, u_y, u_z;
    int id, globe_ind;

    tensor Sigma;
    tensor SigmaTOT;
    SigmaTOT.xx = SigmaTOT.yy = SigmaTOT.zz = 0.0;

for (size_t z = 0 ; z < nz_loc ; z++) {
        for (size_t y = 0 ; y < ny_loc ; y++) {
            for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = xsum + x;
	int y_global = ysum + y;
	int z_global = zsum + z;
				
	size_t ip = ipx(x, y, z,nx_loc,ny_loc);
	globe_ind = (x_global + nx * y_global + nx * ny * z_global);
	id = Id[ip/6];
				
            u_x = p[ip]   - x_global * gridStep;
            u_y = p[ip + 1] - y_global * gridStep;
            u_z = p[ip + 2] - z_global * gridStep;

            if (fabs(u_x) < 1e-10) u_x = 0.0;
            if (fabs(u_y) < 1e-10) u_y = 0.0;
            if (fabs(u_z) < 1e-10) u_z = 0.0;

            computeStress(ipx(x,y,z,nx_loc,ny_loc)/6,x,y,z, Sigma);
                
            if (fabs(Sigma.xx) < 1e-10) Sigma.xx = 0.0;
            if (fabs(Sigma.yy) < 1e-10) Sigma.yy = 0.0;
            if (fabs(Sigma.zz) < 1e-10) Sigma.zz = 0.0;
				
	    if (fabs(Sigma.xy) < 1e-10) Sigma.xy = 0.0;
            if (fabs(Sigma.xz) < 1e-10) Sigma.xz = 0.0;
            if (fabs(Sigma.yz) < 1e-10) Sigma.yz = 0.0;
				
	    if (fabs(Sigma.yx) < 1e-10) Sigma.yx = 0.0;
            if (fabs(Sigma.zx) < 1e-10) Sigma.zx = 0.0;
            if (fabs(Sigma.zy) < 1e-10) Sigma.zy = 0.0;

	    	SigmaTOT.xx += Sigma.xx;
		SigmaTOT.yy += Sigma.yy;
		SigmaTOT.zz += Sigma.zz;
				
            fprintf(sortie, "%i %i %i %i %i %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n", globe_ind,id,x_global,y_global,z_global,u_x,u_y,u_z,Sigma.xx, Sigma.yy, Sigma.zz,Sigma.xy,Sigma.yx,Sigma.xz,Sigma.zx,Sigma.yz,Sigma.zy);

            }
        }
    }

	fclose(sortie);
}



void initialize_rho(double ri, int idp, int nx_loc, int ny_loc, int nz_loc){

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == idp) rho[ipx(x,y,z,nx_loc,ny_loc)/6] = ri;
	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] != idp) rho[ipx(x,y,z,nx_loc,ny_loc)/6] = 0.0;

	}	
	}
	} 

	MPI_Sendrecv (&(rho[0]),1,EW_faceETA,NeighBor[W],flag,&(rho[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[1 * (nx_loc - 1)]),1,EW_faceETA,NeighBor[E],flag,&(rho[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(rho[0]),1,bktypeETA,NeighBor[F],flag,&(rho[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[0]),1,frtypeETA,NeighBor[B],flag,&(rho[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[0]),1,tptypeETA,NeighBor[N],flag,&(rho[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[0]),1,bttypeETA,NeighBor[S],flag,&(rho[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);

}


void DFT_ObjFunc(double tol, double Temp, double mu, int nx_loc, int ny_loc, int nz_loc,int nx_fill, int ny_fill, int nz_fill, int rrkk, int xls, int xrs, int yls, int yrs, int zds, int zus){
	
	double GP, err_tot, rho_tot;
	double err = 1000.0;
	double errsum_loc = 0.0;
	double den = nx * ny * nz;
  	 bool res_cond;

    	double* rhom;
	rhom = new double[nx_fill * ny_fill * nz_fill];

	while (err > tol){

	for (size_t i = 0 ; i < (nx_fill * ny_fill * nz_fill) ; i++) rhom[i] = rho[i];

	double rho_loc = 0.0;
	for (size_t i = 0 ; i < (nx_loc * ny_loc * nz_loc) ; i++) rho_loc += rhom[i];

 	for (size_t z = 0 ; z < nz_loc ; z++){					
	for (size_t y = 0 ; y < ny_loc ; y++){
	for (size_t x = 0 ; x < nx_loc ; x++){

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);

          res_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) && ( z_global >= zds ) && ( z_global <  (nz) - zus) ? false: true; 

	GP = 0.0;
	
	for (int l = 0 ; l < nodes[ipx(x, y, z,nx_loc,ny_loc)/6].nblg ; ++l) GP += rhom[nodes[ipx(x, y, z,nx_loc,ny_loc)/6].ilg[l]/6] + yvar * (1.0 - Eta[nodes[ipx(x, y, z,nx_loc,ny_loc)/6].ilg[l]/6]);
	for (int l = HLG ; l < HLG + nodes[ipx(x, y, z,nx_loc,ny_loc)/6].nbolg ; ++l) GP += rhom[nodes[ipx(x, y, z,nx_loc,ny_loc)/6].ilg[l]/6] + yvar * (1.0 - Eta[nodes[ipx(x, y, z,nx_loc,ny_loc)/6].ilg[l]/6]);
	rho[ipx(x, y, z,nx_loc,ny_loc)/6] = Eta[ipx(x, y, z,nx_loc,ny_loc)/6] / (1.0 + exp( (-1.0/Boltz/Temp) * (mu + wff * GP) ) );
	//if (res_cond) rho[ipx(x,y,z,nx_loc,ny_loc)/6] = 0.999;
	}
	}
	}



	MPI_Sendrecv (&(rho[0]),1,EW_faceETA,NeighBor[W],flag,&(rho[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[1 * (nx_loc - 1)]),1,EW_faceETA,NeighBor[E],flag,&(rho[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(rho[0]),1,bktypeETA,NeighBor[F],flag,&(rho[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[0]),1,frtypeETA,NeighBor[B],flag,&(rho[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[0]),1,tptypeETA,NeighBor[N],flag,&(rho[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[0]),1,bttypeETA,NeighBor[S],flag,&(rho[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);

	for (size_t i = 0 ; i < (nx_loc * ny_loc * nz_loc) ; i++) errsum_loc += pow( (rhom[i] - rho[i]) , 2);
	
	MPI_Barrier(comm3d);
	MPI_Allreduce(&errsum_loc, &err_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&rho_loc, &rho_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	errsum_loc = 0.0;
	err = err_tot/den;
	}
}

void density_gradient(int IDP, int IDS, double trho[], double drx[], double dry[], double drz[], double gs){

	int p0,p1,p2,p3,p4,p5,p6;

	double Dx1, Dy1, Dz1, Dx0, Dy0, Dz0, Dx, Dy, Dz;
	bool xulF, yulF, xllF, yllF, zulF, zllF;
	bool xulT, yulT, xllT, yllT, zulT, zllT;
	bool xug, xlg, yug, ylg, zug, zlg;

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {
		
	int x_global = (xsum + x);
	int y_global = (ysum + y);
	int z_global = (zsum + z);

	if (Id[ipx(x,y,z,nx_loc,ny_loc)/6] == IDP){

	xug = (x_global == nx - 1) ? false : true;
	xlg = (x_global == 0) ? false : true;
	yug = (y_global == ny - 1) ? false : true;
	ylg = (y_global == 0) ? false : true;
	zug = (z_global == nz - 1) ? false : true;
	zlg = (z_global == 0) ? false : true;

	xulF = (x == nx_loc - 1) ? false : true;//if false - continue
	yulF = (y == ny_loc - 1) ? false : true;//if false - continue 
	zulF = (z == nz_loc - 1) ? false : true;//if false - continue
	xllF = (x == 0) ? false : true; //if false - continue 
	yllF = (y == 0) ? false : true; //if false - continue
	zllF = (z == 0) ? false : true; //if false - continue

	xulT = (x == nx_loc - 1) ? true : false;//if true - continue
	yulT = (y == ny_loc - 1) ? true : false;//if true - continue 
	zulT = (z == nz_loc - 1) ? true : false;//if ture - continue
	xllT = (x == 0) ? true : false; //if true - continue 
	yllT = (y == 0) ? true : false; //if true - continue
	zllT = (z == 0) ? true : false; //if true - continue

	if (xulF)	p1 = Id[ipx(x + 1,y,z,nx_loc,ny_loc)/6];
	if (xulT)	p1 = Id[ (nx_loc * ny_loc * nz_loc + z * (ny_loc) + y) ];

	if (yulF)	p2 = Id[ipx(x,y+1,z,nx_loc,ny_loc)/6];
	if (yulT)	p2 = Id[ (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc * 2 + z * (nx_loc + 2) + (x+1)) ];

	if (zulF)	p3 = Id[ipx(x,y,z+1,nx_loc,ny_loc)/6];
	if (zulT)	p3 = Id[(nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * (nx_loc + 2) * nz_loc + (nx_loc + 2) * (ny_loc + 2) + (y + 1) * (nx_loc + 2) + (x + 1))];
		
	if (xllF)	p4 = Id[ipx(x-1,y,z,nx_loc,ny_loc)/6];
	if (xllT)	p4 = Id[(nx_loc * ny_loc * nz_loc + ny_loc * nz_loc + z * ny_loc + (y))];
				
	if (yllF)	p5 = Id[ipx(x, y - 1,z,nx_loc,ny_loc)/6];
	if (yllT)	p5 = Id[(nx_loc * ny_loc * nz_loc + ny_loc * nz_loc * 2 + (nx_loc + 2) * nz_loc + z * (nx_loc + 2) + (x+1))];
				
	if (zllF)	p6 = Id[ipx(x,y,z-1,nx_loc,ny_loc)/6];
	if (zllT)	p6 = Id[(nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * (nx_loc + 2) * nz_loc + (y+1) * (nx_loc + 2) + (x+1))];

	if (xulF)	Dx1 = trho[ipx(x + 1,y,z,nx_loc,ny_loc)/6];
	if (xulT)	Dx1 = trho[ (nx_loc * ny_loc * nz_loc + z * (ny_loc) + y) ];

	if (yulF)	Dy1 = trho[ipx(x,y+1,z,nx_loc,ny_loc)/6];
	if (yulT)	Dy1 = trho[ (nx_loc * ny_loc * nz_loc + ny_loc * nz_loc * 2 + z * (nx_loc + 2) + (x+1)) ];

	if (zulF)	Dz1 = trho[ipx(x,y,z+1,nx_loc,ny_loc)/6];
	if (zulT)	Dz1 = trho[(nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * (nx_loc + 2) * nz_loc + (nx_loc + 2) * (ny_loc + 2) + (y + 1) * (nx_loc + 2) + (x + 1))];
		
	if (xllF)	Dx0 = trho[ipx(x-1,y,z,nx_loc,ny_loc)/6];
	if (xllT)	Dx0 = trho[(nx_loc * ny_loc * nz_loc + ny_loc * nz_loc + z * ny_loc + (y))];
				
	if (yllF)	Dy0 = trho[ipx(x, y - 1,z,nx_loc,ny_loc)/6];
	if (yllT)	Dy0 = trho[(nx_loc * ny_loc * nz_loc + ny_loc * nz_loc * 2 + (nx_loc + 2) * nz_loc + z * (nx_loc + 2) + (x+1))];
				
	if (zllF)	Dz0 = trho[ipx(x,y,z-1,nx_loc,ny_loc)/6];
	if (zllT)	Dz0 = trho[(nx_loc * ny_loc * nz_loc + 2 * ny_loc * nz_loc + 2 * (nx_loc + 2) * nz_loc + (y+1) * (nx_loc + 2) + (x+1))];	

	if (p1 == IDP && p4 == IDP) Dx = (Dx1 - Dx0) / (2.0*gs);
	if (p1 == IDP && p4 == IDS) Dx = (Dx1 - trho[ipx(x,y,z,nx_loc,ny_loc)/6]) / (gs);
	if (p1 == IDS && p4 == IDP) Dx = (trho[ipx(x,y,z,nx_loc,ny_loc)/6] - Dx0) / (gs);

	if (p2 == IDP && p5 == IDP) Dy = (Dy1 - Dy0) / (2.0*gs);
	if (p2 == IDP && p5 == IDS) Dy = (Dy1 - trho[ipx(x,y,z,nx_loc,ny_loc)/6]) / (gs);
	if (p2 == IDS && p5 == IDP) Dy = (trho[ipx(x,y,z,nx_loc,ny_loc)/6] - Dy0) / (gs);

	if (p3 == IDP && p6 == IDP) Dz = (Dz1 - Dz0) / (2.0*gs);
	if (p3 == IDP && p6 == IDS) Dz = (Dz1 - trho[ipx(x,y,z,nx_loc,ny_loc)/6]) / (gs);
	if (p3 == IDS && p6 == IDP) Dz = (trho[ipx(x,y,z,nx_loc,ny_loc)/6] - Dz0) / (gs);

	drx[ipx(x,y,z,nx_loc,ny_loc)/6] = Dx;
	dry[ipx(x,y,z,nx_loc,ny_loc)/6] = Dy;
	drz[ipx(x,y,z,nx_loc,ny_loc)/6] = Dz;

	}

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == IDS){
	drx[ipx(x,y,z,nx_loc,ny_loc)/6] = 0.0;
	dry[ipx(x,y,z,nx_loc,ny_loc)/6] = 0.0;
	drz[ipx(x,y,z,nx_loc,ny_loc)/6] = 0.0;
	}
		
	}
	}
	}
}


void write_output_id(int num, int nx_loc, int ny_loc, int nz_loc){
	
	int globe_ind;
 	FILE * sortie;
 	char nomfic[256];
 	sprintf(nomfic, "id_%i", num);
 	sortie = fopen(nomfic, "w+");

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);
				
	globe_ind = (x_global + nx * y_global + nx * ny * z_global);

	fprintf(sortie,"%i %i\n",globe_ind,Id[ipx(x,y,z,nx_loc,ny_loc)/6]);
	
	}	
	}
	}

fclose(sortie);

}

void write_output_sat(const string& _name, int num, int rrkk, int nx_loc, int ny_loc, int nz_loc){
	
	ostringstream str1,str2;
	str1 << num ; 
	str2 << rrkk ; 
	string app1 = str1.str();
	string app2 = str2.str();
	// string result = _name + to_string(num);
	string result = _name + app1;
	result.append("_");
	// result = result + to_string(rrkk);
	result = result + app2;

	const char * c = result.c_str();
 	FILE * sortie;
 	sortie = fopen(c, "w+");

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);
				
	int globe_ind = (x_global + nx * y_global + nx * ny * z_global);

	fprintf(sortie,"%i %g\n",globe_ind, rho[ipx(x,y,z,nx_loc,ny_loc)/6]);
	
	}	
	}
	}

fclose(sortie);

}

void write_output_results(const string& _name, double a1, double a2, double a3, double a4, double a5, int a6, double a7, double a8){

	const char * c = _name.c_str();
 	FILE * sortie;
 	sortie = fopen(c, "a");

	fprintf(sortie,"%g %g %g %g %g %i %g %g\n",a1,a2,a3,a4,a5,a6,a7,a8);
	fclose(sortie);

}

void write_output_sim_param(double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9, double a10, double a11, double a12, double a13, double a14, double a15, double a16, double a17, double a18, double a19){

FILE * sortie;
sortie = fopen("sim_param", "a");
fprintf(sortie,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19);
fclose(sortie);

}

void write_output_force(int num, int nx_loc, int ny_loc, int nz_loc){
	
	int globe_ind;
 	FILE * sortie;
 	char nomfic[256];
 	sprintf(nomfic, "force-%i", num);
 	sortie = fopen(nomfic, "w+");

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);
				
	globe_ind = (x_global + nx * y_global + nx * ny * z_global);

	fprintf(sortie,"%i %g %g %g\n",globe_ind,force[ipx(x,y,z,nx_loc,ny_loc)],force[ipx(x,y,z,nx_loc,ny_loc) + 1],force[ipx(x,y,z,nx_loc,ny_loc) + 2]);
	
	}	
	}
	}

fclose(sortie);

}

double post_processing(int idp, double tmp[], double mu_sat , int xls, int xrs, int yls, int yrs, int zds, int zus){

	double GP_loc, rhof_loc;
	GP_loc = rhof_loc = 0.0;

	double GP_loc_A, rhof_loc_A, rhof_loc_B;
	GP_loc_A = rhof_loc_A = rhof_loc_B = 0.0;

	double count = 0.0;
          bool res_cond, box_cond;

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = (xsum + x);
	int y_global = (ysum + y);
	int z_global = (zsum + z);

          res_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) && ( z_global >= zds ) && ( z_global <  (nz) - zus) ? true: false;//if true -> continue 
          box_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) && ( z_global >= zds ) && ( z_global <  (nz) - zus) ? false: true;//if false -> continue 

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == idp){

	GP_loc += Boltz * Temp * (rho[ipx(x,y,z,nx_loc,ny_loc)/6] * log(rho[ipx(x,y,z,nx_loc,ny_loc)/6]) + (Eta[ipx(x,y,z,nx_loc,ny_loc)/6] - rho[ipx(x,y,z,nx_loc,ny_loc)/6]) * log(Eta[ipx(x,y,z,nx_loc,ny_loc)/6] - 
	rho[ipx(x,y,z,nx_loc,ny_loc)/6]));
		
	for (int l = 0 ; l < nodes[ipx(x, y, z,nx_loc,ny_loc)/6].nblg ; ++l){
	GP_loc += -wff * rho[nodes[ipx(x, y, z,nx_loc,ny_loc)/6].ilg[l]/6] * rho[ipx(x,y,z,nx_loc,ny_loc)/6];
	}
		
	for (int l = HLG ; l < HLG + nodes[ipx(x, y, z,nx_loc,ny_loc)/6].nbolg ; ++l){
	GP_loc += -wff * rho[nodes[ipx(x, y, z,nx_loc,ny_loc)/6].ilg[l]/6] * rho[ipx(x,y,z,nx_loc,ny_loc)/6];
	}
		
	GP_loc -= (mu - mu_sat) * rho[ipx(x,y,z,nx_loc,ny_loc)/6];
	rhof_loc += rho[ipx(x,y,z,nx_loc,ny_loc)/6];
	
	}

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == idp && res_cond){

	GP_loc_A += Boltz * Temp * (rho[ipx(x,y,z,nx_loc,ny_loc)/6] * log(rho[ipx(x,y,z,nx_loc,ny_loc)/6]) + (Eta[ipx(x,y,z,nx_loc,ny_loc)/6] - rho[ipx(x,y,z,nx_loc,ny_loc)/6]) * log(Eta[ipx(x,y,z,nx_loc,ny_loc)/6] - rho[ipx(x,y,z,nx_loc,ny_loc)/6]));
		
	for (int l = 0 ; l < nodes[ipx(x, y, z,nx_loc,ny_loc)/6].nblg ; ++l){
	GP_loc_A += -wff * rho[nodes[ipx(x, y, z,nx_loc,ny_loc)/6].ilg[l]/6] * rho[ipx(x,y,z,nx_loc,ny_loc)/6];
	}
		
	for (int l = HLG ; l < HLG + nodes[ipx(x, y, z,nx_loc,ny_loc)/6].nbolg ; ++l){
	GP_loc_A += -wff * rho[nodes[ipx(x, y, z,nx_loc,ny_loc)/6].ilg[l]/6] * rho[ipx(x,y,z,nx_loc,ny_loc)/6];
	}
		
	GP_loc_A -= (mu - mu_sat) * rho[ipx(x,y,z,nx_loc,ny_loc)/6];
	rhof_loc_A += rho[ipx(x,y,z,nx_loc,ny_loc)/6];
	count = count + 1.0;
	}

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == idp && box_cond) rhof_loc_B += rho[ipx(x,y,z,nx_loc,ny_loc)/6];

	}
	}
	}

	double GP_tot, rhof_tot;
	double GP_tot_A, rhof_tot_A, rhof_tot_B;
	double count_g, countB_g;
	
	MPI_Barrier(comm3d);
	MPI_Allreduce(&GP_loc, &GP_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&rhof_loc, &rhof_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&GP_loc_A, &GP_tot_A, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&rhof_loc_A, &rhof_tot_A, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&rhof_loc_B, &rhof_tot_B, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&count, &count_g, 1, MPI_DOUBLE, MPI_SUM, comm3d);

	tmp[0] = rhof_tot;//rhof   box+reservoir
	tmp[1] = rhof_tot_A;//rhof   box
	tmp[2] = rhof_tot_B;

	tmp[3] = GP_tot;//local GP box+reservoir
	tmp[4] = GP_tot_A / count_g;//local GP box
  	return GP_tot;
}


void stress_calculation(size_t k, double nn, double MU, double gs, double Temp, tensor & stress, double trho[], double drx[], double dry[], double drz[], double mu_sat ){

	double p0 = 0.0;
	double p1 = 0.0;
	double cn = 6.0;
	
	p0 += -Boltz * Temp * (	trho[k] * log(trho[k]) + (1.0 - trho[k]) * log(1.0 - trho[k])	) + (mu) * trho[k];
	p0 += (cn * wff / 2.0) * pow(	trho[k],2	);

	stress.xx =  p0 - (cn * pow(gs,2) * wff / 4.0	) * (	drx[k] * drx[k] + dry[k] * dry[k] + drz[k] * drz[k]	) + (cn * pow(gs,2) / 2.0 * wff * (	drx[k] * drx[k]	)	);
	stress.xy =  cn * pow(gs,2) / 2.0 * wff * drx[k] * dry[k];
	stress.xz =  cn * pow(gs,2) / 2.0 * wff * drx[k] * drz[k];
	stress.yy =  p0 - (cn * pow(gs,2) * wff / 4.0	) * (	drx[k] * drx[k] + dry[k] * dry[k] + drz[k] * drz[k]	) + (cn * pow(gs,2) / 2.0 * wff * (	dry[k] * dry[k]	)	);
	stress.yz =  cn * pow(gs,2) / 2.0 * wff * dry[k] * drz[k];
	stress.zz =  p0 - (cn * pow(gs,2) * wff / 4.0	) * (	drx[k] * drx[k] + dry[k] * dry[k] + drz[k] * drz[k]	) + (cn * pow(gs,2) / 2.0 * wff * (	drz[k] * drz[k]	)	);

	stress.xx = stress.xx/pow(gs,3);
	stress.xy = stress.xy/pow(gs,3);
	stress.xz = stress.xz/pow(gs,3);
	stress.yy = stress.yy/pow(gs,3);
	stress.yz = stress.yz/pow(gs,3);
	stress.zz = stress.zz/pow(gs,3);

	stress.yx = stress.xy;
	stress.zx = stress.xz;
	stress.zy = stress.yz;

}

void write_output_fluid_stress_data(const string& _name, int num, int rrkk, double NearNeigh, double ChemPot, double LatSpac, double Temp, int idp, int ids, double trho[], double drx[], double dry[], double drz[], double mu_sat, int xls, int xrs, int yls, int yrs, int zds, int zus){
	
	tensor sigma, sFres;
	bool res_cond; 
	double count, sigxx,sigyy,sigzz,sigxy,sigxz,sigyz;
	sigxx = sigyy = sigzz = sigxy = sigxz = sigyz = count = 0.0;

	ostringstream str1,str2;
	str1 << num ; 
	str2 << rrkk ; 
	string app1 = str1.str();
	string app2 = str2.str();
	// string result = _name + to_string(num);
	string result = _name + app1;
	result.append("_");
	// result = result + to_string(rrkk);
	result = result + app2;

	const char * c = result.c_str();
 	FILE * sortie;
 	sortie = fopen(c, "w+");

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = (xsum + x);
	int y_global = (ysum + y);
	int z_global = (zsum + z);

          res_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) && ( z_global >= zds ) && ( z_global <  (nz) - zus) ? false: true;//if false -> continue 
				
	int globe_ind = (x_global + nx * y_global + nx * ny * z_global);
	stress_calculation(ipx(x,y,z,nx_loc,ny_loc)/6, NearNeigh, ChemPot, LatSpac, Temp, sigma, trho, drx, dry, drz, mu_sat);

	if(res_cond){
	count++;
	sigxx += sigma.xx;
	sigyy += sigma.yy;
	sigzz += sigma.zz;
	sigxy += sigma.xy;
	sigxz += sigma.xz;
	sigyz += sigma.yz;
	}

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == idp) fprintf(sortie,"%i %.6g\n",globe_ind,(sigma.xx+sigma.yy+sigma.zz)/3.0);
	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == ids) fprintf(sortie,"%i %.6g\n",globe_ind,0.0);
			
	}	
	}
	}

	fclose(sortie);

	double count_tot,sigxx_tot,sigyy_tot,sigzz_tot,sigxy_tot,sigxz_tot,sigyz_tot;
	MPI_Barrier(comm3d);
	MPI_Allreduce(&count, &count_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&sigxx, &sigxx_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&sigyy, &sigyy_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&sigzz, &sigzz_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&sigxy, &sigxy_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&sigxz, &sigxz_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	MPI_Allreduce(&sigyz, &sigyz_tot, 1, MPI_DOUBLE, MPI_SUM, comm3d);
	sFres.xx = sigxx_tot / count_tot;
	sFres.yy = sigyy_tot / count_tot;
	sFres.zz = sigzz_tot / count_tot;
	sFres.xy = sigxy_tot / count_tot;
	sFres.xz = sigxz_tot / count_tot;
	sFres.yz = sigyz_tot / count_tot;

	if (rrkk == 0){
	FILE * sortieA;
	sortieA = fopen("sigF_res","a");
	fprintf(sortieA,"%g %g %g %g %g %g\n",sFres.xx,sFres.yy,sFres.zz,sFres.xy,sFres.xz,sFres.yz);
	fclose(sortieA);
	}

}

/*
void fluid_stress_computation(double NearNeigh, double ChemPot, double LatSpac, double Temp, int idp, int ids, double trho[], double drx[], double dry[], double drz[], double mu_sat ){
	
	tensor sigma;
	int globe_ind;

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = (xsum + x);
	int y_global = (ysum + y);
	int z_global = (zsum + z);
				
	globe_ind = (x_global + nx * y_global + nx * ny * z_global);

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == idp){
	stress_calculation(ipx(x,y,z,nx_loc,ny_loc)/6, NearNeigh, ChemPot, LatSpac, Temp, sigma, trho, drx, dry, drz, mu_sat);

	sigF[ipx(x,y,z,nx_loc,ny_loc)] = sigma.xx;
	sigF[ipx(x,y,z,nx_loc,ny_loc) + 1] = sigma.xy;
	sigF[ipx(x,y,z,nx_loc,ny_loc) + 2] = sigma.xz;
	sigF[ipx(x,y,z,nx_loc,ny_loc) + 3] = sigma.yy;
	sigF[ipx(x,y,z,nx_loc,ny_loc) + 4] = sigma.yz;
	sigF[ipx(x,y,z,nx_loc,ny_loc) + 5] = sigma.zz;
		
	}
					
	}	
	}
	}
}
*/

/*
void write_output_press(const string& _name, int num, int rrkk, double NearNeigh, double ChemPot, double LatSpac, double Temp, int idp, int ids, double trho[], double drx[], double dry[], double drz[], double mu_sat){
	
	tensor sigma;
	double* stress_decomp;
	stress_decomp = new double[2];

	string result = _name + to_string(num);
	result.append("_");
	result = result + to_string(rrkk);
	const char * c = result.c_str();
 	FILE * sortie;
 	sortie = fopen(c, "w+");

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = (xsum + x);
	int y_global = (ysum + y);
	int z_global = (zsum + z);
				
	int globe_ind = (x_global + nx * y_global + nx * ny * z_global);
	stress_calculation(ipx(x,y,z,nx_loc,ny_loc)/6, NearNeigh, ChemPot, LatSpac, Temp, sigma, trho, drx, dry, drz, mu_sat, stress_decomp);

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == idp) fprintf(sortie,"%i %i %.6g\n",globe_ind,Id[ipx(x,y,z,nx_loc,ny_loc)/6],(sigma.xx+sigma.yy+sigma.zz)/3.0);
	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == ids) fprintf(sortie,"%i %i %.6g\n",globe_ind,Id[ipx(x,y,z,nx_loc,ny_loc)/6],0.0);
				
	}	
	}
	}

	fclose(sortie);

	delete [] stress_decomp;

}*/

/*
void write_output_stress(const string& _name, int num, int rrkk, double NearNeigh, double ChemPot, double LatSpac, double Temp, int idp, int ids, int nx_loc, int ny_loc, int nz_loc){
	
	tensor sigma;
	int globe_ind;

	string result = _name + to_string(num);
	result.append("_");
	result = result + to_string(rrkk);
	const char * c = result.c_str();
 	FILE * sortie;
 	sortie = fopen(c, "w+");

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);
				
	globe_ind = (x_global + nx * y_global + nx * ny * z_global);

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == idp){
	stress_calculation(ipx(x,y,z,nx_loc,ny_loc)/6, NearNeigh, ChemPot, LatSpac, Temp, sigma);
	fprintf(sortie,"%i %i %.6g %.6g %.6g %.6g %.6g %.6g\n",globe_ind,Id[ipx(x,y,z,nx_loc,ny_loc)/6],sigma.xx,sigma.xy,sigma.xz,sigma.yy,sigma.yz,sigma.zz);
	sigF[ipx(x,y,z,nx_loc,ny_loc)] = sigma.xx;
	sigF[ipx(x,y,z,nx_loc,ny_loc) + 1] = sigma.xy;
	sigF[ipx(x,y,z,nx_loc,ny_loc) + 2] = sigma.xz;
	sigF[ipx(x,y,z,nx_loc,ny_loc) + 3] = sigma.yy;
	sigF[ipx(x,y,z,nx_loc,ny_loc) + 4] = sigma.yz;
	sigF[ipx(x,y,z,nx_loc,ny_loc) + 5] = sigma.zz;		
	}
		
	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == ids){
	fprintf(sortie,"%i %i %.6g %.6g %.6g %.6g %.6g %.6g\n",globe_ind,Id[ipx(x,y,z,nx_loc,ny_loc)/6],0.0,0.0,0.0,0.0,0.0,0.0);
	}
			
	}	
	}
	}

fclose(sortie);

}
*/

void write_output_data(const string& _name, int num, int rrkk, double data[] ){
	
	ostringstream str1,str2;
	str1 << num ; 
	str2 << rrkk ; 
	string app1 = str1.str();
	string app2 = str2.str();
	// string result = _name + to_string(num);
	string result = _name + app1;
	result.append("_");
	// result = result + to_string(rrkk);
	result = result + app2;

	const char * c = result.c_str();
 	FILE * sortie;
 	sortie = fopen(c, "w+");

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = (xsum + x);
	int y_global = (ysum + y);
	int z_global = (zsum + z);

	int globe_ind = (x_global + nx * y_global + nx * ny * z_global);

	fprintf( sortie,"%i %g\n",globe_ind, data[ipx(x,y,z,nx_loc,ny_loc)/6]);
			
	}	
	}
	}

fclose(sortie);

}

void ApplyFluidForce(int idp, int ids, int rrkk, double NearNeigh, double ChemPot, double LatSpac, double Temp){
	
	tensor sigma;
	double fx,fy,fz, fx_tot, fy_tot, fz_tot, count;
	int tag;

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	count = fx_tot = fy_tot = fz_tot = fx = fy = fz = 0.0;
	
	tag = 0;

	if(Id[ipx(x,y,z,nx_loc,ny_loc)/6] == ids){
	
	for (int l = 0 ; l < nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nb ; ++l) {
	if(Id[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l]/6] == idp){
	tag = tag + 1;
	count = count + 1.0;

	fx = sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l]] * (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) + 
					   sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 1] * (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) +
				           sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 2] * (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
	fy = sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 1] * (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) + 
					   sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 3] * (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) +
				           sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 4] * (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
	fz = sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 2] * (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) + 
					   sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 4] * (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) +
				           sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 5] * (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);

	fx_tot += fx;
	fy_tot += fy;
	fz_tot += fz;

	force[ipx(x,y,z,nx_loc,ny_loc)] += fx;
	force[ipx(x,y,z,nx_loc,ny_loc) + 1] += fy;
	force[ipx(x,y,z,nx_loc,ny_loc) + 2] += fz;

	}
	}

	for (int l = HALF ; l < HALF + nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nbo ; ++l) {
	if(Id[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l]/6] == idp){
	tag = tag + 1;
	count = count + 1.0;

	fx = sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l]] * (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) + 
					   sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 1] * (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) +
				           sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 2] * (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
	fy = sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 1] * (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) + 
					   sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 3] * (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) +
				           sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 4] * (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);
	fz = sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 2] * (_nx[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) + 
					   sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 4] * (_ny[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]) +
				           sigF[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l] + 5] * (_nz[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]] / _L0[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].I[l]]);

	fx_tot += fx;
	fy_tot += fy;
	fz_tot += fz;

	force[ipx(x,y,z,nx_loc,ny_loc)] += fx;
	force[ipx(x,y,z,nx_loc,ny_loc) + 1] += fy;
	force[ipx(x,y,z,nx_loc,ny_loc) + 2] += fz;

	}
	}

	}

	if(tag > 0){

	force[ipx(x,y,z,nx_loc,ny_loc)] = fx_tot/count;
	force[ipx(x,y,z,nx_loc,ny_loc) + 1] = fy_tot/count;
	force[ipx(x,y,z,nx_loc,ny_loc) + 2] = fz_tot/count;

	}

	}	
	}
	}
}


void write_output_sig_solid(int num){

	tensor Sigma;

 	FILE * sortie;
 	char nomfic[256];
 	sprintf(nomfic, "sigS-%i", num);
 	sortie = fopen(nomfic, "w+");
	int globe_ind;
	
	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);
				
	size_t ip = ipx(x, y, z,nx_loc,ny_loc);
	globe_ind = (x_global + nx * y_global + nx * ny * z_global);
        computeStress(ipx(x,y,z,nx_loc,ny_loc)/6,x,y,z, Sigma);
	if (fabs(Sigma.xx) < 1e-10) Sigma.xx = 0.0;
        if (fabs(Sigma.yy) < 1e-10) Sigma.yy = 0.0;
        if (fabs(Sigma.zz) < 1e-10) Sigma.zz = 0.0;
			
	if (fabs(Sigma.xy) < 1e-10) Sigma.xy = 0.0;
        if (fabs(Sigma.xz) < 1e-10) Sigma.xz = 0.0;
        if (fabs(Sigma.yz) < 1e-10) Sigma.yz = 0.0;
				
	if (fabs(Sigma.xy) < 1e-10) Sigma.yx = 0.0;
        if (fabs(Sigma.xz) < 1e-10) Sigma.zx = 0.0;
        if (fabs(Sigma.yz) < 1e-10) Sigma.zy = 0.0;
	
	fprintf(sortie,"%i %g %g %g %g %g %g\n",globe_ind,Sigma.xx,Sigma.xy,Sigma.xz,Sigma.yy,Sigma.yz,Sigma.zz);

	}	
	}
	}

fclose(sortie);

}

void write_output_disp(int num){

	double u_x, u_y, u_z;
 	FILE * sortie;
 	char nomfic[256];
 	sprintf(nomfic, "disp-%i", num);
 	sortie = fopen(nomfic, "w+");
	int globe_ind;
	
	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);
				
	size_t ip = ipx(x, y, z,nx_loc,ny_loc);
	globe_ind = (x_global + nx * y_global + nx * ny * z_global);

	u_x = p[ipx(x,y,z,nx_loc,ny_loc)] - x_global * gridStep;
	u_y = p[ipx(x,y,z,nx_loc,ny_loc) + 1] - y_global * gridStep;
	u_z = p[ipx(x,y,z,nx_loc,ny_loc) + 2] - z_global * gridStep;

	fprintf(sortie,"%i %g  %g  %g\n",globe_ind,u_x,u_y,u_z);

	}	
	}
	}

fclose(sortie);

}

void ApplyPressureOut(int idp, int ids, double Force, int rrkk){

 	for (size_t z = 0 ; z < nz_loc ; z++){					
		for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {
				
				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);
										
		if (Id[ipx(x,y,z,nx_loc,ny_loc)/6] == ids){
			
		for (int l = 0 ; l < nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nb ; ++l) {
		if (Id[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l]/6] == idp){

		force[ipx(x,y,z,nx_loc,ny_loc)] += Force;
		force[ipx(x,y,z,nx_loc,ny_loc) + 1] += Force;
		force[ipx(x,y,z,nx_loc,ny_loc) + 2] += Force;
		}
		}
		for (int l = HALF ; l < HALF + nodes[ipx(x,y,z,nx_loc,ny_loc)/6].nbo ; ++l) {
		if (Id[nodes[ipx(x,y,z,nx_loc,ny_loc)/6].ipx[l]/6] == idp){

		force[ipx(x,y,z,nx_loc,ny_loc)] += Force;
		force[ipx(x,y,z,nx_loc,ny_loc) + 1] += Force;
		force[ipx(x,y,z,nx_loc,ny_loc) + 2] += Force;	
		}
		}
		}

		}
		}
	}
}

void  reservoir_creation(int idincl, double ideta, int xls, int xrs, int yls, int yrs, int zds, int zus){

   bool res_cond;

   for (size_t z = 0 ; z < nz_loc ; z++) {					
   for (size_t y = 0 ; y < ny_loc ; y++) {
   for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);
		
          res_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) && ( z_global >= zds ) && ( z_global <  (nz) - zus) ? false: true;//if true -> continue 

   if (res_cond){
   Id[ipx(x,y,z,nx_loc,ny_loc)/6] = idincl;
   Eta[ipx(x,y,z,nx_loc,ny_loc)/6] = ideta;
   }

   }
   }
   }

	MPI_Barrier(comm3d);
	MPI_Sendrecv (&(Id[0]),1,EW_faceID,NeighBor[W],flag,&(Id[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_INT,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[1 * (nx_loc - 1)]),1,EW_faceID,NeighBor[E],flag,&(Id[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_INT,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(Id[0]),1,bktypeID,NeighBor[F],flag,&(Id[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_INT,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,frtypeID,NeighBor[B],flag,&(Id[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_INT,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,tptypeID,NeighBor[N],flag,&(Id[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_INT,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(Id[0]),1,bttypeID,NeighBor[S],flag,&(Id[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_INT,NeighBor[N],flag,comm3d,&status);

	MPI_Sendrecv (&(Eta[0]),1,EW_faceETA,NeighBor[W],flag,&(Eta[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[1 * (nx_loc - 1)]),1,EW_faceETA,NeighBor[E],flag,&(Eta[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(Eta[0]),1,bktypeETA,NeighBor[F],flag,&(Eta[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,frtypeETA,NeighBor[B],flag,&(Eta[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,tptypeETA,NeighBor[N],flag,&(Eta[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(Eta[0]),1,bttypeETA,NeighBor[S],flag,&(Eta[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);

}

void  initialize_reservoir(double ideta, int xls, int xrs, int yls, int yrs, int zds, int zus){

   bool res_cond;

   for (size_t z = 0 ; z < nz_loc ; z++) {					
   for (size_t y = 0 ; y < ny_loc ; y++) {
   for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = (xsum + x);
				int y_global = (ysum + y);
				int z_global = (zsum + z);

          res_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) && ( z_global >= zds ) && ( z_global <  (nz) - zus) ? false: true;//if true -> continue 
  
 if (res_cond){
   rho[ipx(x,y,z,nx_loc,ny_loc)/6] = 1.0;
   }

   }
   }
   }
	MPI_Sendrecv (&(rho[0]),1,EW_faceETA,NeighBor[W],flag,&(rho[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[1 * (nx_loc - 1)]),1,EW_faceETA,NeighBor[E],flag,&(rho[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(rho[0]),1,bktypeETA,NeighBor[F],flag,&(rho[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[0]),1,frtypeETA,NeighBor[B],flag,&(rho[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[0]),1,tptypeETA,NeighBor[N],flag,&(rho[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(rho[0]),1,bttypeETA,NeighBor[S],flag,&(rho[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);

}


void correct_rho(double tnrho[], double lb, double ub){

	for (size_t z = 0 ; z < nz_loc ; z++) {
	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {
	double lrho = rho[ipx(x,y,z,nx_loc,ny_loc)/6];
	tnrho[ipx(x,y,z,nx_loc,ny_loc)/6] = lrho;
	if(lrho < lb) tnrho[ipx(x,y,z,nx_loc,ny_loc)/6] = 1e-20;
	if(lrho > ub) tnrho[ipx(x,y,z,nx_loc,ny_loc)/6] = 0.9999999;

	}	
	}
	} 

	MPI_Sendrecv (&(tnrho[0]),1,EW_faceETA,NeighBor[W],flag,&(tnrho[nz_loc*nx_loc*ny_loc*1]), (nz_loc * ny_loc * 1), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(tnrho[1 * (nx_loc - 1)]),1,EW_faceETA,NeighBor[E],flag,&(tnrho[nz_loc*nx_loc*ny_loc*1 + nz_loc * ny_loc * 1]), nz_loc * ny_loc * 1, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(tnrho[0]),1,bktypeETA,NeighBor[F],flag,&(tnrho[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(tnrho[0]),1,frtypeETA,NeighBor[B],flag,&(tnrho[nz_loc * nx_loc * ny_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1]), ((nx_loc + 2) * nz_loc * 1), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(tnrho[0]),1,tptypeETA,NeighBor[N],flag,&(tnrho[nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * nz_loc * 1 * 2]),(nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(tnrho[0]),1,bttypeETA,NeighBor[S],flag,&(tnrho[nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + 2 * (nx_loc + 2) * nz_loc * 1 + nx_fill * ny_fill * 1]), (nx_fill * ny_fill * 1),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);

}


int main(int argc, char *argv[]){

	readSimParam();

	int nx_lt_shift = Lx_lt_shift / gridStep + 1;
	int nx_rt_shift = Lx_rt_shift / gridStep + 1;
	int ny_lt_shift = Ly_lt_shift / gridStep + 1;
	int ny_rt_shift = Ly_rt_shift / gridStep + 1;
	int nz_up_shift = Lz_up_shift / gridStep + 1;
	int nz_dn_shift = Lz_dn_shift / gridStep + 1;

	if (Lx_lt_shift == 0) nx_lt_shift = 0;
	if (Lx_rt_shift == 0) nx_rt_shift = 0;
	if (Ly_lt_shift == 0) ny_lt_shift = 0;
	if (Ly_rt_shift == 0) ny_rt_shift = 0;
	if (Lz_up_shift == 0) nz_up_shift = 0;
	if (Lz_dn_shift == 0) nz_dn_shift = 0;

	/*nx = Lnx / gridStep + 1 + nx_lt_shift + nx_rt_shift;
	ny = Lny / gridStep + 1 + ny_lt_shift + ny_rt_shift;
	nz = Lnz / gridStep + 1 + nz_dn_shift + nz_up_shift;*/

	/*int nx_lt_shift = Lx_lt_shift / gridStep + 1;
	int nx_rt_shift = Lx_rt_shift / gridStep + 1;
	int ny_lt_shift = 0;
	int ny_rt_shift = 0;
	int nz_dn_shift = 0;
	int nz_up_shift = 0;*/
	
	nx = Lnx / gridStep + 1 + nx_lt_shift + nx_rt_shift;
	ny = Lny / gridStep + 1;
	nz = Lnz / gridStep + 1;

	MPI_Comm comm;
	int iter;
	double fret;
	int nproc, rank;     
    	int reorder = 0;
	
	int ndims = 3;
	ndof = 6;
	int dims[ndims];
	int periods[ndims];
	int coord[ndims];

	MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm,&nproc);
    MPI_Comm_rank(comm,&rank);

	lengthX = new int[nxp];
	lengthY = new int[nyp];
	lengthZ = new int[nzp];
		
	decompose1d(lengthX,nx,nxp);
	decompose1d(lengthY,ny,nyp);
	decompose1d(lengthZ,nz,nzp);
	
	periods[0] = 1;
    	periods[1] = 1;
	periods[2] = 1;
	
    dims[0] = nxp;
    dims[1] = nyp;
	dims[2] = nzp;

    MPI_Cart_create(comm, ndims, dims, periods, reorder, &comm3d);
	MPI_Cart_get(comm3d , ndims , dims , periods , coord);

	NeighBor[0] = MPI_PROC_NULL;
    NeighBor[1] = MPI_PROC_NULL;
    NeighBor[2] = MPI_PROC_NULL;
    NeighBor[3] = MPI_PROC_NULL;
	NeighBor[4] = MPI_PROC_NULL;
    NeighBor[5] = MPI_PROC_NULL;

    MPI_Cart_shift(comm3d,0,1,&NeighBor[W],&NeighBor[E]);
	MPI_Cart_shift(comm3d,1,1,&NeighBor[F],&NeighBor[B]);
  	MPI_Cart_shift(comm3d,2,1,&NeighBor[S],&NeighBor[N]);
	
	nx_loc = lengthX[coord[0]];
	ny_loc = lengthY[coord[1]];
	nz_loc = lengthZ[coord[2]];
		
	nx_fill = (nx_loc + 2);
	ny_fill = (ny_loc + 2);
	nz_fill = (nz_loc + 2);
	
	xsum = 0;
	ysum = 0;
	zsum = 0;
	int xx,yy,zz;
	
	for(xx = 0 ; xx < coord[0] ; xx++) xsum += lengthX[xx];
	for(yy = 0 ; yy < coord[1] ; yy++) ysum += lengthY[yy];
	for(zz = 0 ; zz < coord[2] ; zz++) zsum += lengthZ[zz];
			
	p = new double[nx_fill * ny_fill * nz_fill * ndof];
	sigF = new double[nx_fill * ny_fill * nz_fill * ndof];
	
	for(size_t i = 0 ; i < (nx_fill * ny_fill * nz_fill * ndof) ; i++){
	p[i] = 0.0;
	sigF[i] = 0.0;
	}
	
	int count = 0;
	for (int k = 0 ; k < nz_loc ; k++) {
	for (int j = 0 ; j < ny_loc ; j++) {
	for (int i = 0 ; i < nx_loc ; i++) {
	int x = xsum + i;
	int y = ysum + j;
	int z = zsum + k;
	p[ipx(i,j,k,nx_loc,ny_loc)] = x * gridStep;
	p[ipx(i,j,k,nx_loc,ny_loc) + 1] = y * gridStep;
	p[ipx(i,j,k,nx_loc,ny_loc) + 2] = z * gridStep;
	count++;
	}
	}
	}
	
	force = new double[nx_fill * ny_fill * nz_fill * ndof];
	
	for(size_t i = 0 ; i < (nx_fill * ny_fill * nz_fill * ndof) ; i++) force[i] = 0.0;
	
	Id = new int[nx_fill * ny_fill * nz_fill];
	Eta = new double[nx_fill * ny_fill * nz_fill];
	
	blocklen = new int[3 * nz_loc];
	int ip = 0;
	for(size_t i = 0 ; i < (nz_loc) ; i++){
	blocklen[ip++] = ndof;
	blocklen[ip++] = ndof * nx_loc;
	blocklen[ip++] = ndof;
	}
	
	blocklenID = new int[3 * nz_loc];
	ip = 0;
	for(size_t i = 0 ; i < (nz_loc) ; i++){
	blocklenID[ip++] = 1;
	blocklenID[ip++] = 1 * nx_loc;
	blocklenID[ip++] = 1;
	}
			
	disp_b = new int[3 * nz_loc];
	ip = 0;
	for(size_t i = 0 ; i < (nz_loc) ; i++){
	disp_b[ip++] = nx_loc * ny_loc * nz_loc * ndof + ny_loc * nz_loc * ndof + i * (ny_loc) * ndof;
	disp_b[ip++] = i * nx_loc * ny_loc * ndof;
	disp_b[ip++] = nx_loc * ny_loc * nz_loc * ndof + i * (ny_loc) * ndof;
	}
	
	disp_bID = new int[3 * nz_loc];
	ip = 0;
	for(size_t i = 0 ; i < (nz_loc) ; i++){
	disp_bID[ip++] = nx_loc * ny_loc * nz_loc * 1 + ny_loc * nz_loc * 1 + i * (ny_loc) * 1;
	disp_bID[ip++] = i * nx_loc * ny_loc * 1;
	disp_bID[ip++] = nx_loc * ny_loc * nz_loc * 1 + i * (ny_loc) * 1;
	}
			
	disp_f = new int[3 * nz_loc];
	ip = 0;
	for(size_t i = 0 ; i < (nz_loc) ; i++){
	disp_f[ip++] = nx_loc * ny_loc * nz_loc * ndof + ny_loc * nz_loc * ndof + (ny_loc-1) * ndof + i * ny_loc * ndof;
	disp_f[ip++] = nx_loc * (ny_loc - 1) * ndof + i * nx_loc * ny_loc * ndof;
	disp_f[ip++] = nx_loc * ny_loc * nz_loc * ndof + (ny_loc - 1) * ndof + i * ny_loc * ndof;
	}
	
	disp_fID = new int[3 * nz_loc];
	ip = 0;
	for(size_t i = 0 ; i < (nz_loc) ; i++){
	disp_fID[ip++] = nx_loc * ny_loc * nz_loc * 1 + ny_loc * nz_loc * 1 + (ny_loc-1) * 1 + i * ny_loc * 1;
	disp_fID[ip++] = nx_loc * (ny_loc - 1) * 1 + i * nx_loc * ny_loc * 1;
	disp_fID[ip++] = nx_loc * ny_loc * nz_loc * 1 + (ny_loc - 1) * 1 + i * ny_loc * 1;
	}
			
	blocklen_dt = new int[3 * ny_loc + 2];
	blocklen_dt [0] = (nx_loc + 2) * ndof;
	blocklen_dt [3 * ny_loc + 1] = (nx_loc + 2) * ndof;
	ip = 1;
	for(size_t i = 0 ; i < (ny_loc) ; i++){
	blocklen_dt[ip++] = ndof;
	blocklen_dt[ip++] = ndof * nx_loc;
	blocklen_dt[ip++] = ndof;
	}
	
	blocklen_dtID = new int[3 * ny_loc + 2];
	blocklen_dtID [0] = (nx_loc + 2) * 1;
	blocklen_dtID [3 * ny_loc + 1] = (nx_loc + 2) * 1;
	ip = 1;
	for(size_t i = 0 ; i < (ny_loc) ; i++){
	blocklen_dtID[ip++] = 1;
	blocklen_dtID[ip++] = 1 * nx_loc;
	blocklen_dtID[ip++] = 1;
	}
		
	disp_tp = new int[3 * ny_loc + 2];
	disp_tp [3 * ny_loc + 1] = nx_loc * ny_loc * nz_loc * ndof + 2 * ny_loc * nz_loc * ndof + (nx_loc + 2) * (nz_loc - 1) * ndof;
	disp_tp [0] = nx_loc * ny_loc * nz_loc * ndof + 2 * nz_loc * ny_loc * ndof + (nx_loc + 2) * (nz_loc) * ndof + (nx_loc + 2) * (nz_loc - 1) * ndof;
	ip = 1;
	for(size_t i = 0 ; i < (ny_loc) ; i++){
	disp_tp[ip++] = nx_loc * ny_loc * nz_loc * ndof + ny_loc * nz_loc * ndof + ny_loc * (nz_loc - 1) * ndof + i * ndof;
	disp_tp[ip++] = nx_loc * ny_loc * (nz_loc - 1) * ndof + i * nx_loc * ndof;
	disp_tp[ip++] = nx_loc * ny_loc * nz_loc * ndof + ny_loc * (nz_loc - 1) * ndof + i * ndof;
	}
	
	disp_tpID = new int[3 * ny_loc + 2];
	disp_tpID [3 * ny_loc + 1] = nx_loc * ny_loc * nz_loc * 1 + 2 * ny_loc * nz_loc * 1 + (nx_loc + 2) * (nz_loc - 1) * 1;
	disp_tpID [0] = nx_loc * ny_loc * nz_loc * 1 + 2 * nz_loc * ny_loc * 1 + (nx_loc + 2) * (nz_loc) * 1 + (nx_loc + 2) * (nz_loc - 1) * 1;
	ip = 1;
	for(size_t i = 0 ; i < (ny_loc) ; i++){
	disp_tpID[ip++] = nx_loc * ny_loc * nz_loc * 1 + ny_loc * nz_loc * 1 + ny_loc * (nz_loc - 1) * 1 + i * 1;
	disp_tpID[ip++] = nx_loc * ny_loc * (nz_loc - 1) * 1 + i * nx_loc * 1;
	disp_tpID[ip++] = nx_loc * ny_loc * nz_loc * 1 + ny_loc * (nz_loc - 1) * 1 + i * 1;
	}
		
	disp_bt = new int[3 * ny_loc + 2];
	disp_bt [0] = nx_loc * ny_loc * nz_loc * ndof + ny_loc * nz_loc * ndof * 2 + (nx_loc + 2) * nz_loc * ndof;
	disp_bt [3 * ny_loc + 1] = nx_loc * ny_loc * nz_loc * ndof + ny_loc * nz_loc * ndof * 2;
	ip = 1;
	for(size_t i = 0 ; i < (ny_loc) ; i++){
	disp_bt[ip++] = nx_loc * ny_loc * nz_loc * ndof + ny_loc * nz_loc * ndof + i * ndof;
	disp_bt[ip++] = i * nx_loc * ndof;
	disp_bt[ip++] = nx_loc * ny_loc * nz_loc * ndof + i * ndof;
	}
	
	disp_btID = new int[3 * ny_loc + 2];
	disp_btID [0] = nx_loc * ny_loc * nz_loc * 1 + ny_loc * nz_loc * 1 * 2 + (nx_loc + 2) * nz_loc * 1;
	disp_btID [3 * ny_loc + 1] = nx_loc * ny_loc * nz_loc * 1 + ny_loc * nz_loc * 1 * 2;
	ip = 1;
	for(size_t i = 0 ; i < (ny_loc) ; i++){
	disp_btID[ip++] = nx_loc * ny_loc * nz_loc * 1 + ny_loc * nz_loc * 1 + i * 1;
	disp_btID[ip++] = i * nx_loc * 1;
	disp_btID[ip++] = nx_loc * ny_loc * nz_loc * 1 + i * 1;
	}
	
	MPI_Type_vector(nz_loc * ny_loc, ndof, ndof * nx_loc, MPI_DOUBLE, &EW_face);  
	MPI_Type_commit(&EW_face);
	
	MPI_Type_vector(nz_loc * ny_loc, 1, 1 * nx_loc, MPI_INT, &EW_faceID);  
	MPI_Type_commit(&EW_faceID);
	
	MPI_Type_vector(nz_loc * ny_loc, 1, 1 * nx_loc, MPI_DOUBLE, &EW_faceETA);  
	MPI_Type_commit(&EW_faceETA);
	
	MPI_Type_vector(nz_loc, ndof * nx_loc, nx_loc * ny_loc * ndof, MPI_DOUBLE, &BF_face);  
	MPI_Type_commit(&BF_face);
	
	MPI_Type_vector(nz_loc, 1 * nx_loc, nx_loc * ny_loc * 1, MPI_INT, &BF_faceID);  
	MPI_Type_commit(&BF_faceID);
	
	MPI_Type_vector(nz_loc, 1 * nx_loc, nx_loc * ny_loc * 1, MPI_DOUBLE, &BF_faceETA);  
	MPI_Type_commit(&BF_faceETA);

	MPI_Type_indexed(3*nz_loc,blocklen,disp_b,MPI_DOUBLE,&bktype);
	MPI_Type_commit(&bktype);
	
	MPI_Type_indexed(3*nz_loc,blocklenID,disp_bID,MPI_INT,&bktypeID);
	MPI_Type_commit(&bktypeID);
	
	MPI_Type_indexed(3*nz_loc,blocklenID,disp_bID,MPI_DOUBLE,&bktypeETA);
	MPI_Type_commit(&bktypeETA);
	
	MPI_Type_indexed(3*nz_loc,blocklen,disp_f,MPI_DOUBLE,&frtype);
	MPI_Type_commit(&frtype);
	
	MPI_Type_indexed(3*nz_loc,blocklenID,disp_fID,MPI_INT,&frtypeID);
	MPI_Type_commit(&frtypeID);
	
	MPI_Type_indexed(3*nz_loc,blocklenID,disp_fID,MPI_DOUBLE,&frtypeETA);
	MPI_Type_commit(&frtypeETA);
		
	MPI_Type_indexed(3 * ny_loc + 2,blocklen_dt,disp_tp,MPI_DOUBLE,&tptype);
	MPI_Type_commit(&tptype);
	
	MPI_Type_indexed(3 * ny_loc + 2,blocklen_dtID,disp_tpID,MPI_INT,&tptypeID);
	MPI_Type_commit(&tptypeID);
	
	MPI_Type_indexed(3 * ny_loc + 2,blocklen_dtID,disp_tpID,MPI_DOUBLE,&tptypeETA);
	MPI_Type_commit(&tptypeETA);	

	MPI_Type_indexed(3*ny_loc+2,blocklen_dt,disp_bt,MPI_DOUBLE,&bttype);
	MPI_Type_commit(&bttype);
	
	MPI_Type_indexed(3*ny_loc+2,blocklen_dtID,disp_btID,MPI_INT,&bttypeID);
	MPI_Type_commit(&bttypeID);
	
	MPI_Type_indexed(3*ny_loc+2,blocklen_dtID,disp_btID,MPI_DOUBLE,&bttypeETA);
	MPI_Type_commit(&bttypeETA);

	int IDP = 0;
	int IDS = 1;
	
	_coef = new double[Q];
	_coefB = new double[Q];
	_coefT = new double[Q];

	save_MPI_info(rank,nx_loc,ny_loc,nz_loc,coord[0],coord[1],coord[2],xsum,lengthX[coord[0]],ysum,lengthY[coord[1]],zsum,lengthZ[coord[2]]);	
		
	init(gridStep, ndof, nx_loc, ny_loc, nz_loc, xsum, ysum, zsum);
	for(size_t i = 0 ; i < (nx_fill * ny_fill * nz_fill) ; i++){
	Id[i] = IDS;	
	Eta[i] = 0.0;
	}

	rho = new double[nx_fill * ny_fill * nz_fill];
	double* drhoX;
	double* drhoY;
	double* drhoZ;
	double* nrho;

	drhoX = new double[nx_loc * ny_loc * nz_loc];
	drhoY = new double[nx_loc * ny_loc * nz_loc];
	drhoZ = new double[nx_loc * ny_loc * nz_loc];
	nrho = new double[nx_fill * ny_fill * nz_fill];
	
	for (size_t i = 0 ; i < (nx_loc * ny_loc * nz_loc) ; i++){
	drhoX[i] = 0.0;
	drhoY[i] = 0.0;
	drhoZ[i] = 0.0;
	nrho [i] = 0.0;
	}

	//Import_Structure("FinalCoor.txt", IDS, 0.0, nu_part,nx_loc, ny_loc, nz_loc, xsum, ysum, zsum, Lx_lt_shift, Ly_lt_shift, Lz_dn_shift);
	cylinder(IDP, 1.0, XC, YC, ZC, Rpore, nx_loc, ny_loc, nz_loc, xsum, ysum, zsum, nx_lt_shift, ny_lt_shift, nz_dn_shift);
	reservoir_creation(IDP,1.0, nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift,nz_dn_shift,nz_up_shift);
	write_output_id(rank,nx_loc,ny_loc,nz_loc);
	
	double NN = 6.0;
	double surf_ten = 0.072;// N*m^-1
	Boltz = 1.38064852 * pow(10,-23);//m^2*kg*s^-2*K^-1
    	//Boltz = 8.3144621;//J*mol^-1*K^-1
	wff = Boltz * Tcw * 4.0 / NN ;//J*mol^-1//m^2*kg*s^-2=N.m=J
	wmf = yvar * wff;//J*mol^-1//m^2*kg*s^-2=N.m=J
  	Temp = Tstar * wff / Boltz;
	
	double mu_sat = -wff * NN / 2.0;
	double vp0 = exp(mu_sat/Boltz/Temp);
	double a_lg = sqrt(wff/2.0/surf_ten);//m

	if (rank == 0) write_output_sim_param(Boltz,wff,wmf,yvar,Tcw,Tstar,Temp,mu_sat,vp0,gridStep,Lnx,Lny,Lnz, Lx_lt_shift, Lx_rt_shift, Ly_lt_shift, Ly_rt_shift, Lz_dn_shift, Lz_up_shift);	

	//non-dimensional values
	Temp = Boltz*Temp/wff;
	wmf = wmf/wff;
	Boltz = Boltz/Boltz;
	mu_sat = mu_sat / wff;
	wff = wff / wff;
	
	double h,vp;
	double* res;
	double* NoPore;
	res = new double[5];
	NoPore = new double[3];
		
	h = 0.001;
	mu = mu_sat + Boltz * Temp * log(h);
	double rhoInit = exp((1/Boltz/Temp) * mu);
	initialize_rho(0.001,IDP,nx_loc,ny_loc,nz_loc);
    	//initialize_reservoir(0.999,nx_lt_shift,nx_rt_shilsft,ny_lt_shift,ny_rt_shift,nz_dn_shift,nz_up_shift);

	count = 0;
	for (int i = 1; i < 101; ++i){
	double step = static_cast<double>(i) / 100.0;
	h = step;
	mu = mu_sat + Boltz * Temp * log(h);
	// mu = (-8.0 + step) * wff;
	// h = exp( (mu - mu_sat) / Boltz / Temp);
	vp = exp(mu/Boltz/Temp);
	mu = mu / wff;
	
	DFT_ObjFunc(obj_tol, Temp, mu, nx_loc, ny_loc, nz_loc,nx_fill,ny_fill,nz_fill,rank,nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift,nz_dn_shift,nz_up_shift);	
	double GrandPot = post_processing(IDP,res, 0, nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift,nz_dn_shift,nz_up_shift);
	double phi = vol_frac(NoPore,IDP, rank, nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift,nz_dn_shift,nz_up_shift);
	
	correct_rho(nrho,rho_thresh_lb, rho_thresh_ub);
	density_gradient(IDP,IDS,nrho,drhoX,drhoY,drhoZ,gridStep);

	if (rank == 0) write_output_results("sim_data_tot",res[0]/NoPore[0],res[1]/NoPore[1],res[2]/NoPore[2],mu,h,vp/vp0,res[3]/NoPore[0],res[4]);
	if (rank == 0) cout<<"step = " <<count<<" mu: "<<mu<<endl;
		
	write_output_fluid_stress_data("f_",count, rank, NN, mu, gridStep,  Temp, IDP, IDS, nrho, drhoX, drhoY, drhoZ,mu_sat, nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift,nz_dn_shift,nz_up_shift);
	write_output_data("sat_",count,rank,rho);
	count++;
	}

	/*
	double NNodes = nx - 1;
	double Mind = 4000.0;
	double nu = 0.2499;
	Elastic_Input_ISO(Mind, nu, NNodes);
	
	Ktable.resize(3);for (size_t i = 0 ; i < 3 ; i++) Ktable[i].resize(3);
	
	Ktable[0][0] = 0.0;
	Ktable[1][1] = Mind;
	Ktable[2][2] = 0.0;
	Ktable[0][1] = 0.0;
	Ktable[1][0] = 0.0;
	Ktable[2][0] = 0.0;
	Ktable[0][2] = 0.0;
	Ktable[2][1] = 0.0;
	Ktable[1][2] = 0.0;
		
	MPI_Sendrecv (&(p[0]),1,EW_face,NeighBor[W],flag,&(p[nz_loc*nx_loc*ny_loc*ndof]), (nz_loc * ny_loc * ndof), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(p[ndof * (nx_loc - 1)]),1,EW_face,NeighBor[E],flag,&(p[nz_loc*nx_loc*ny_loc*ndof + nz_loc * ny_loc * ndof]), nz_loc * ny_loc * ndof, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(p[0]),1,bktype,NeighBor[F],flag,&(p[nz_loc * nx_loc * ny_loc * ndof + 2 * nz_loc * ny_loc * ndof]), ((nx_loc + 2) * nz_loc * ndof), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(p[0]),1,frtype,NeighBor[B],flag,&(p[nz_loc * nx_loc * ny_loc * ndof + 2 * nz_loc * ny_loc * ndof + (nx_loc + 2) * nz_loc * ndof]), ((nx_loc + 2) * nz_loc * ndof), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(p[0]),1,tptype,NeighBor[N],flag,&(p[nx_loc * ny_loc * nz_loc * ndof + 2 * nz_loc * ny_loc * ndof + (nx_loc + 2) * nz_loc * ndof * 2]),(nx_fill * ny_fill * ndof),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(p[0]),1,bttype,NeighBor[S],flag,&(p[nx_loc * ny_loc * nz_loc * ndof + 2 * ny_loc * nz_loc * ndof + 2 * (nx_loc + 2) * nz_loc * ndof + nx_fill * ny_fill * ndof]), (nx_fill * ny_fill * ndof),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);

	MPI_Sendrecv (&(sigF[0]),1,EW_face,NeighBor[W],flag,&(sigF[nz_loc*nx_loc*ny_loc*ndof]), (nz_loc * ny_loc * ndof), MPI_DOUBLE,NeighBor[E],flag,comm3d,&status);
	MPI_Sendrecv (&(sigF[ndof * (nx_loc - 1)]),1,EW_face,NeighBor[E],flag,&(sigF[nz_loc*nx_loc*ny_loc*ndof + nz_loc * ny_loc * ndof]), nz_loc * ny_loc * ndof, MPI_DOUBLE,NeighBor[W],flag,comm3d,&status);	
	MPI_Sendrecv (&(sigF[0]),1,bktype,NeighBor[F],flag,&(sigF[nz_loc * nx_loc * ny_loc * ndof + 2 * nz_loc * ny_loc * ndof]), ((nx_loc + 2) * nz_loc * ndof), MPI_DOUBLE,NeighBor[B],flag,comm3d,&status);
	MPI_Sendrecv (&(sigF[0]),1,frtype,NeighBor[B],flag,&(sigF[nz_loc * nx_loc * ny_loc * ndof + 2 * nz_loc * ny_loc * ndof + (nx_loc + 2) * nz_loc * ndof]), ((nx_loc + 2) * nz_loc * ndof), MPI_DOUBLE,NeighBor[F],flag,comm3d,&status);
	MPI_Sendrecv (&(sigF[0]),1,tptype,NeighBor[N],flag,&(sigF[nx_loc * ny_loc * nz_loc * ndof + 2 * nz_loc * ny_loc * ndof + (nx_loc + 2) * nz_loc * ndof * 2]),(nx_fill * ny_fill * ndof),MPI_DOUBLE,NeighBor[S],flag,comm3d,&status);
	MPI_Sendrecv (&(sigF[0]),1,bttype,NeighBor[S],flag,&(sigF[nx_loc * ny_loc * nz_loc * ndof + 2 * ny_loc * nz_loc * ndof + 2 * (nx_loc + 2) * nz_loc * ndof + nx_fill * ny_fill * ndof]), (nx_fill * ny_fill * ndof),MPI_DOUBLE,NeighBor[N],flag,comm3d,&status);
	
	force_init();
 	//FabricTensor(IDPa, rank);
 	//ApplyPorePressure(IDPa, IDS, res[5]/NoPore[1], NoPore[1], rank);
	ApplyFluidForce(IDPa,IDS,rank,NN,mu,a_lg,Temp);
	ApplyPressureOut(IDPb,IDS, 0.0, rank);
	Force_BC(IDPa,IDS,rank,nx_loc, ny_loc, nz_loc, xsum, ysum, zsum);
	Force_BC(IDPb,IDS,rank,nx_loc, ny_loc, nz_loc, xsum, ysum, zsum);
	frprmn(p, (nx_fill * ny_fill * nz_fill * ndof), ftol, 1500, &iter, &fret, Up_springs, Grad_Up_springs_app,rank);

	write_output_force(rank,nx_loc,ny_loc,nz_loc);
	write_output_sig_solid(rank);
	write_output_disp(rank);
	*/
	
delete [] Id;
delete [] p;
delete [] sigF;
delete [] force;
delete [] blocklen;
delete [] disp_b;
delete [] disp_f;
delete [] blocklen_dt;
delete [] disp_tp;
delete [] disp_bt;
delete [] blocklenID;
delete [] disp_bID;
delete [] disp_fID;
delete [] blocklen_dtID;
delete [] disp_tpID;
delete [] disp_btID;
delete [] lengthX;
delete [] lengthY;
delete [] lengthZ;
delete [] _coef;
delete [] _coefB;
delete [] _coefT;
delete [] Eta;
delete [] rho;
delete [] drhoX;
delete [] drhoY;
delete [] drhoZ;
delete [] res;
delete [] NoPore;
delete [] nrho;

MPI_Type_free(&EW_face);
MPI_Type_free(&EW_faceID);
MPI_Type_free(&EW_faceETA);
	
MPI_Type_free(&BF_face);
MPI_Type_free(&BF_faceID);
MPI_Type_free(&BF_faceETA);
	
MPI_Type_free(&bktype);
MPI_Type_free(&bktypeID);
MPI_Type_free(&bktypeETA);
	
MPI_Type_free(&frtype);
MPI_Type_free(&frtypeID);
MPI_Type_free(&frtypeETA);
	
MPI_Type_free(&tptype);
MPI_Type_free(&tptypeID);
MPI_Type_free(&tptypeETA);
	
MPI_Type_free(&bttype);
MPI_Type_free(&bttypeID);
MPI_Type_free(&bttypeETA);
	
MPI_Comm_free(&comm3d);

MPI_Finalize();

return 0;
}
