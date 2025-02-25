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

using namespace std;

int nx, ny, nz, nxp, nyp, nzp, nx_og, ny_og, nz_og, Lnx, Lny, Lnz;
double gridStep, type, strain, nit, ftol, yvar, XC, YC, ZC, Rpore, Tcw, Tstar;
double Lx_lt_shift, Lx_rt_shift, Ly_lt_shift, Ly_rt_shift, Lz_up_shift, Lz_dn_shift;
double rho_thresh_lb, rho_thresh_ub, obj_tol;
size_t nu_part;


int ipx(int x, int y, int z, int xsize, int ysize){
	return (x + xsize * y + xsize * ysize * z);
}

void  Import_MPI_info(const string& _name, int tmp[]){
	
	
	int rrkk, nx_loc, ny_loc, nz_loc, xcoor, ycoor, zcoor, xs, xsp, ys, ysp, zs, zsp;
	fstream file(_name);

	file >> rrkk >> nx_loc >> ny_loc >> nz_loc >> xcoor >> ycoor >> zcoor >> xs >> xsp >> ys >> ysp >> zs >> zsp;
	tmp[0] = rrkk;
	tmp[1] = nx_loc;
	tmp[2] = ny_loc;
	tmp[3] = nz_loc;
	tmp[4] = xcoor;
	tmp[5] = ycoor;
	tmp[6] = zcoor;
	tmp[7] = xs;
	tmp[8] = xsp;
	tmp[9] = ys;
	tmp[10] = ysp;
	tmp[11] = zs;
	tmp[12] = zsp;
}


void  Import_sat(const string& _name, int num, int rrkk, int nxl, int nyl, int nzl, double tmp[]){
	
	ostringstream str1,str2;
	str1 << num ; 
	str2 << rrkk ; 
	string app1 = str1.str();
	string app2 = str2.str();
	// string result = _name + to_string(num);
	string resultIN = _name + app1;
	resultIN.append("_");
	// result = result + to_string(rrkk);
	resultIN = resultIN + app2;
	
	int globe_ind;
	double rho;
	
	fstream file(resultIN);

	for (size_t k = 0 ; k < (nxl * nyl * nzl) ; k++){
	file >> globe_ind >> rho;
	tmp[globe_ind] = rho;
	}
}

void  Import_id(const string& _name, int rrkk, int nxl, int nyl, int nzl, int tmp[]){
	
	ostringstream str1;
	str1 << rrkk;
	string app1 = str1.str();
	string resultIN = _name + app1;
	
	int globe_ind, tid;
	
	fstream file(resultIN);

	for (size_t k = 0 ; k < (nxl * nyl * nzl) ; k++){
	file >> globe_ind >> tid;
	tmp[globe_ind] = tid;
	}
}


void  Export(const string& _name, int num, double tmp[], size_t tot){
	
	ostringstream str1;
	str1 << num;
	string app1 = str1.str();
	string resultIN = _name + app1;
	
   const char * c = resultIN.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");	
   for (size_t k = 0 ; k < (tot) ; k++) fprintf(sortie,"%g\n",tmp[k]);
   fclose(sortie);
}

void  Export_id(const string& _name, int tmp[], size_t tot){
	
   const char * c = _name.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");	

   for (size_t k = 0 ; k < (tot) ; k++) fprintf(sortie,"%i\n",tmp[k]);
   fclose(sortie);
}

void  frame_change_id(int tmpA[], int tmpB[], int xls, int xrs, int yls, int yrs, int zds, int zus){

   int countA = 0;
   int countB = 0;
   bool res_cond;

   for (int z = 0 ; z < nz ; z++) {
   for (int y = 0 ; y < ny ; y++) {
   for (int x = 0 ; x < nx ; x++) {

          res_cond = ( x >= xls ) && ( x <  (nx) - xrs ) && ( y >= yls ) && ( y < (ny) - yrs ) && ( z >= zds ) && ( z <  (nz) - zus) ? true: false;//if true -> continue 
   
   if (res_cond){
   tmpB[countB] = tmpA[countA];
   countB++;
   }
   countA++;
   }
   }
   }
}

void  frame_change(double tmpA[], double tmpB[], int xls, int xrs, int yls, int yrs, int zds, int zus){

   int countA = 0;
   int countB = 0;
   bool res_cond;

   for (int z = 0 ; z < nz ; z++) {
   for (int y = 0 ; y < ny ; y++) {
   for (int x = 0 ; x < nx ; x++) {

          res_cond = ( x >= xls ) && ( x <  (nx) - xrs ) && ( y >= yls ) && ( y < (ny) - yrs ) && ( z >= zds ) && ( z <  (nz) - zus) ? true: false;//if true -> continue 

   if(res_cond){
   tmpB[countB] = tmpA[countA];
   countB++;
   }
   countA++;
   }
   }
   }
}

void  Import_stress(const string& _name, int num, int rrkk, int nxl, int nyl, int nzl, double tmpXX[], double tmpYY[], double tmpZZ[]){
	
	ostringstream str1,str2;
	str1 << num ; 
	str2 << rrkk ; 
	string app1 = str1.str();
	string app2 = str2.str();
	// string result = _name + to_string(num);
	string resultIN = _name + app1;
	resultIN.append("_");
	// result = result + to_string(rrkk);
	resultIN = resultIN + app2;
	
	int globe_ind, tid;
	double txx,txy,txz,tyy,tyz,tzz;
	
	fstream file(resultIN);

	for (size_t k = 0 ; k < (nxl * nyl * nzl) ; k++){
	file >> globe_ind >> tid >> txx >> txy >> txz >> tyy >> tyz >> tzz;
	tmpXX[globe_ind] = txx;
	//tmpXY[globe_ind] = txy;
	//tmpXZ[globe_ind] = txz;
	tmpYY[globe_ind] = tyy;
	//tmpYZ[globe_ind] = tyz;
	tmpZZ[globe_ind] = tzz;
	}
}

void  Import_press(const string& _name, int num, int rrkk, int nxl, int nyl, int nzl, double tpress[] ){
	
	ostringstream str1,str2;
	str1 << num ; 
	str2 << rrkk ; 
	string app1 = str1.str();
	string app2 = str2.str();
	// string result = _name + to_string(num);
	string resultIN = _name + app1;
	resultIN.append("_");
	// result = result + to_string(rrkk);
	resultIN = resultIN + app2;
	
	int globe_ind, tid;
	double press;
	
	fstream file(resultIN);

	for (size_t k = 0 ; k < (nxl * nyl * nzl) ; k++){
	file >> globe_ind >> tid >> press;
	tpress[globe_ind] = press;
	}
}

void  Import_fluid_data(const string& _name, int num, int rrkk, int nxl, int nyl, int nzl, double tpress[]){
	
	ostringstream str1,str2;
	str1 << num ; 
	str2 << rrkk ; 
	string app1 = str1.str();
	string app2 = str2.str();
	// string result = _name + to_string(num);
	string resultIN = _name + app1;
	resultIN.append("_");
	// result = result + to_string(rrkk);
	resultIN = resultIN + app2;

	int globe_ind, tid;
	double press;
	double sxx, sxx_pa, sxx_pb;
	
	fstream file(resultIN);

	for (size_t k = 0 ; k < (nxl * nyl * nzl) ; k++){
	file >> globe_ind >> press;
	tpress[globe_ind ] = press;
	}
}


void compute_cumulants(const string& _name, double tmp[], int tmpID[], int tot, int flag) 
{
	
	string nameA = _name;
	nameA.append("_c1");
	string nameB = _name;
	nameB.append("_c2");
	string nameC = _name;
	nameC.append("_c3");
	string nameD = _name;
	nameD.append("_c4");
	
	const char * cA = nameA.c_str();
	const char * cB = nameB.c_str();
	const char * cC = nameC.c_str();
	const char * cD = nameD.c_str();

	double m1, m2, count, ExprA, ExprB, ExprC;
	m1 = m2 = count = ExprA = ExprB = ExprC = 0.0;
	
	if (flag == 2){
	for (int i = 0 ; i < tot ; i++){
	count++;
	m1 += tmp[i];
	m2 += pow(tmp[i],2);
	}
	}
	
	if (flag == 1){
	for (int i = 0 ; i < tot ; i++){
	if (tmpID[i] > 0){
	count++;
	m1 += tmp[i];
	m2 += pow(tmp[i],2);		
	}
	}
	}
	
	if (flag == 0){
	for (int i = 0 ; i < tot ; i++){
	if (tmpID[i] == 0){
	count++;
	m1 += tmp[i];
	m2 += pow(tmp[i],2);		
	}
	}
	}

	m1 = m1/(count);
	m2 = m2/(count);
	
	for (int i = 0 ; i < tot ; i++){
	ExprA += pow((tmp[i] - m1),2);
	ExprB += pow((tmp[i] - m1),3);
	ExprC += pow((tmp[i] - m1),4);
	}

	double Mean = m1;//first cumulant
	double Std = sqrt(ExprA/(count - 1));
	//double Var = m2 - pow(m1,2);//second cumulant
	double Var = pow(Std,2);
	double Skew = (ExprB/count) / pow(sqrt(ExprA/count),3);
	double Kurt = (ExprC/count) / pow(ExprA/count,2);
	
	FILE * sortie;
	sortie = fopen(cA,"a");
	fprintf(sortie,"%g\n",Mean);
	fclose(sortie);
	
	FILE * sortieB;
	sortieB = fopen(cB,"a");
	fprintf(sortieB,"%g\n",Var);
	fclose(sortieB);
	
	FILE * sortieC;
	sortieC = fopen(cC,"a");
	fprintf(sortieC,"%g\n",Skew);
	fclose(sortieC);
	
	FILE * sortieD;
	sortieD = fopen(cD,"a");
	fprintf(sortieD,"%g\n",Kurt);
	fclose(sortieD);

}

void histogram(const string& _name, int num, double tmp[], int tot, double val_min, double val_max, double binwidth){
	
	ostringstream str1;
	str1 << num;
	string app1 = str1.str();
	string resultIN = _name + app1;

	// string nameA = _name;
	// nameA.append("_hist_");
	// nameA = nameA + to_string(num);

	const char * c = resultIN.c_str();
	
	if(val_min >= 0.0 && val_max >= 0.0){
	
	int nbins;
	int* hist;
	double* intv;
	nbins = (int)ceil ( (val_max) / binwidth);
	hist = new int[nbins];
	intv = new double[nbins];
	double intv0 = val_min + (binwidth)/2.0;
	for(int i = 0 ; i < nbins ; i++) hist[i] = 0;
	for(int i = 0 ; i < nbins ; i++){
	intv0 += binwidth;
	intv[i] = intv0;
	}

	for (int i = 0 ; i < tot ; i++){
    int bucket = (int)floor(tmp[i] / binwidth);
    hist[bucket - (int)floor(val_min/binwidth) ] += 1;
	}	

	FILE * sortie;
	sortie = fopen(c,"w+");
	for(int i = 0 ; i < nbins ; i++) fprintf(sortie,"%g %i\n",intv[i],hist[i]);
	delete [] hist;
	delete [] intv;

	}
	
	if(val_min < 0.0 && val_max < 0.0){
		
	int nbins;
	int* hist;
	double* intv;
	nbins = (int)ceil ( (-val_min) / binwidth);
	hist = new int[nbins];
	intv = new double[nbins];
	double intv0 = val_min + (binwidth)/2.0;
	for(int i = 0 ; i < nbins ; i++) hist[i] = 0;	
	for(int i = 0 ; i < nbins ; i++){
	intv0 += binwidth;
	intv[i] = intv0;
	}

	for (int i = 0 ; i < tot ; i++){
    int bucket = (int)floor(-tmp[i] / binwidth);
    hist[bucket - (int)floor(-val_max/binwidth)] += 1;
	}
	FILE * sortie;
	sortie = fopen(c,"w+");
	for(int i = 0 ; i < nbins ; i++) fprintf(sortie,"%g %i\n",intv[i],hist[nbins - 1 - i]);
	delete [] hist;
	delete [] intv;

	}
	
	if(val_min < 0.0 && val_max >= 0.0){

	int nbins_pos, nbins_neg;
	int* hist_pos; 
	int* hist_neg;
	double* intv;
   	nbins_pos = (int)ceil( (val_max)/binwidth);
   	nbins_neg = (int)ceil( (-val_min)/binwidth);
		
   	hist_pos = new int[nbins_pos];	
   	hist_neg = new int[nbins_neg];
	intv = new double[nbins_pos + nbins_neg];
   	for(int i = 0 ; i < nbins_pos ; i++) hist_pos[i] = 0;
   	for(int i = 0 ; i < nbins_neg ; i++) hist_neg[i] = 0;
	double intv0 = val_min + (binwidth)/2.0;
	for(int i = 0 ; i < (nbins_pos + nbins_neg) ; i++){
	intv0 += binwidth;
	intv[i] = intv0;
	}

	for (int i = 0 ; i < tot ; i++){
	if(tmp[i] >= 0){
    int bucket = (int)floor(tmp[i] / binwidth);
    hist_pos[bucket] += 1;
	}
	if(tmp[i] < 0){
    int bucket = (int)floor(-tmp[i] / binwidth);
    hist_neg[bucket] += 1;
	}
	}
	FILE * sortie;
	sortie = fopen(c,"w+");
	for(int i = 0 ; i < nbins_neg ; i++) fprintf(sortie,"%g %i\n",intv[i],hist_neg[nbins_neg - 1 - i]);
	for(int i = 0 ; i < nbins_pos ; i++) fprintf(sortie,"%g %i\n",intv[i + nbins_neg],hist_pos[i]);
	delete [] hist_pos;
	delete [] hist_neg;
	delete [] intv;	
	}
}


int volfrac_calculation(int tid[], int nxo, int nyo, int nzo) 
{
	double tot = (nxo - 1) * (nyo - 1) * (nzo - 1);
	int NP = 0;
	int NS = 0;

	for (int z = 0 ; z < (nxo * nyo * nzo) ; z++) { 
	if(tid[z] == 0){NP++;};
	if(tid[z] > 0){NS++;};
	}
	double nump = NP;
	double fp = nump/tot;
	double fs = 1.0 - fp;

	cout<<"fp: "<<fp<<" fs: "<<fs<<endl;
	return NP;
}

double find_min(const string& _name,double tmp[], int tid[], int tot, int flag){

	string nameA = _name;
	nameA.append("_min");
	const char * cA = nameA.c_str();
	
	double small = tmp[0];
	
	if(flag == 2){
	for (int i = 0 ; i < tot ; i++){
	if(tmp[i] < small) small = tmp[i];
	}
	}
	
	if(flag == 1){
	for (int i = 0 ; i < tot ; i++){
	if(tid[i] > 0 && tmp[i] < small) small = tmp[i];
	}
	}
	
	if(flag == 0){
	for (int i = 0 ; i < tot ; i++){
	if(tid[i] == 0 && tmp[i] < small) small = tmp[i];
	}
	}

	FILE * sortie;
	sortie = fopen(cA,"a");
	fprintf(sortie,"%g\n",small);
	fclose(sortie);
	
	return small;
}

double find_max(const string& _name,double tmp[], int tid[], int tot, int flag){

	string nameA = _name;
	nameA.append("_max");
	const char * cA = nameA.c_str();
	
	double large = tmp[0];
	
	if(flag == 2){
	for (int i = 0 ; i < tot ; i++){
	if(tmp[i] > large) large = tmp[i];
	}
	}
	
	if(flag == 1){
	for (int i = 0 ; i < tot ; i++){
	if(tid[i] > 0 && tmp[i] > large) large = tmp[i];
	}
	}
	
	if(flag == 0){
	for (int i = 0 ; i < tot ; i++){
	if(tid[i] == 0 && tmp[i] > large) large = tmp[i];
	}
	}
	
	FILE * sortie;
	sortie = fopen(cA,"a");
	fprintf(sortie,"%g\n",large);
	fclose(sortie);

	return large;
}

void  Import_data(const string& _name, double tdata[], int N){
	
	double m,v,s,k;

	fstream file(_name);
	
	for (int i = 0 ; i < N ; i++){					
	file >> m >> v >> s >> k;
	tdata[i] = m;
 	}
}

/*
double Create_structure(int id[],double tmpX[], double tmpY[], double tmpZ[],double rad[], size_t N){

	double xc, yc, zc, crit, r;	
	double NP = 0.0;
	for (size_t i = 0 ; i < N ; i++) { 

	xc = tmpX[i];
	yc = tmpY[i];
	zc = tmpZ[i];
	r = rad[i];

	for (size_t z = 0 ; z < nz ; z++){
   	for (size_t y = 0 ; y < ny ; y++){
	for (size_t x = 0 ; x < nx ; x++){

    	crit=(1.*x-1.*xc)*(1.*x-1.*xc)+(1.*y-1.*yc)*(1.*y-1.*yc)+(1.*z-1.*zc)*(1.*z-1.*zc);
	if(crit<=1.*r*r){	
	size_t ik = x + nx * y + nx * ny * z;
	id[ik]  = 1; 
	NP++;
	}
   	}
    	}
    	}
	}
	double tot = (nx-1) * (ny-1) * (nz-1);
	double fp = NP/ tot;
	return fp;

}

void  d_dist(double tmpX[], double tmpY[], double tmpZ[], double rad[], size_t N, double D[]){
	
	FILE * sortie;
	sortie = fopen("d_dist","w+");
	double xi,yi,zi,xj,yj,zj;
	size_t count = 0;
	for (size_t i = 0 ; i < N-1 ; i++){	
	xi = tmpX[i];
	yi = tmpY[i];
	zi = tmpZ[i];
	for (size_t j = i+1 ; j < N ; j++){			
	xj = tmpX[j];
	yj = tmpY[j];
	zj = tmpZ[j];
	
	fprintf(sortie,"%g\n",sqrt( pow(xi-xj,2) + pow(yi-yj,2) + pow(zi-zj,2) ));
	D[count++] = sqrt( pow(xi-xj,2) + pow(yi-yj,2) + pow(zi-zj,2));

 	}
	}
	fclose(sortie);
}
*/

void data_allocation(double data_og[], int tid[], int tot, double data_new[], int phase_id){
	int count = 0;
	for (int i = 0 ; i < tot ; i++){
	if(tid[i] == phase_id){
	data_new[count] = data_og[i];
	count++;
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

void readSimData(const string& _name, int tmp[], int TOT){
	
	int a1;
	double a2,a3,a4,a5,a6;
	//const char * _cname = _name.c_str();
	ifstream file(_name);
	for (int k = 0 ; k < TOT ; ++k){
	file >> a2 >> a3 >> a4 >> a5 >> a6 >> a1;
	tmp[k] = a1;
	}
}

void save_vtk_ascii(const string& _name, int id[], int nx_i, int ny_i, int nz_i)
{

    size_t nx_box = (nx_i) ;
    size_t ny_box = (ny_i) ;
    size_t nz_box = (nz_i) ;

   const char * c = _name.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");

    /*FILE * sortie;
    char nomfic[256];
    sprintf(nomfic, "result%04i.vtk", 5);
    sortie = fopen(nomfic, "w");*/

    fprintf(sortie, "# vtk DataFile Version 2.0\n");
    fprintf(sortie, "Sortie domaine LB+LINK\n");
    fprintf(sortie, "ASCII\n");
    fprintf(sortie, "DATASET RECTILINEAR_GRID\n");
    fprintf(sortie, "DIMENSIONS %zu %zu %zu\n", nx_box, ny_box, nz_box);
    
    fprintf(sortie, "X_COORDINATES %zu float\n", nx_box);
    for (size_t i = 0; i <  nx_box; i++) {
        fprintf(sortie, "%.4e ", i * gridStep);
    }
    fprintf(sortie, "\n");
    
    fprintf(sortie, "Y_COORDINATES %zu float\n", ny_box);
    for (size_t i = 0; i < ny_box; i++) {
        fprintf(sortie, "%.4e ", i * gridStep);
    }
    fprintf(sortie, "\n");
    
    fprintf(sortie, "Z_COORDINATES %zu float\n", nz_box);
    for (size_t i = 0; i < nz_box ; i++) {
        fprintf(sortie, "%.4e ", i * gridStep);
    }
    fprintf(sortie, "\n");
    
    fprintf(sortie, "POINT_DATA %zu\n", nx_box * ny_box * nz_box);

    size_t ik = 0;
    fprintf(sortie, "SCALARS MatterId int 1\n");
    fprintf(sortie, "LOOKUP_TABLE default\n");
    for (size_t i = 0 ; i < (nx_box * ny_box * nz_box) ; i++){
    if(id[i] == 0)  fprintf(sortie, "%i\n", 0);
    if(id[i] > 0)  fprintf(sortie, "%i\n", 1);
    }
    
    fprintf(sortie, "VECTORS Displacement float\n");
    for (size_t i = 0 ; i < (nx_box * ny_box * nz_box) ; i++) fprintf(sortie, "%.4e %.4e %.4e\n", 0.0, 0.0, 0.0);
    
    
    fprintf(sortie, "VECTORS StressDiag float\n");
    for (size_t i = 0 ; i < (nx_box * ny_box * nz_box) ; i++) fprintf(sortie, "%.4e %.4e %.4e\n", 0.0, 0.0, 0.0);
    
    fclose(sortie);
}


int main (){
	
   readSimParam();

	nx_og  = Lnx / gridStep + 1;
	ny_og =  Lny / gridStep + 1;
	nz_og =  Lnz / gridStep + 1;

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





   int nprocessors = nxp * nyp * nzp;
   int IDP = 0;
   int IDS = 1;
   double val_max, val_min, binsize;
   int NP,NS, simsample;
	
   char nomfic[256];
   int* locsize;
   locsize = new int[nprocessors * 3];
   int* out;
   out = new int[13];	
   int* rank;
   rank = new int[nprocessors];
   int* xll;
   xll = new int[nprocessors];	
   int* xul;
   xul = new int[nprocessors];	
   int* yll;
   yll = new int[nprocessors];		
   int* yul;
   yul = new int[nprocessors];		
   int* zll;
   zll = new int[nprocessors];	
   int* zul;
   zul = new int[nprocessors];

   int* id;
   int* id_cut;

   double* data_mod;	
   double* sat;
   double* press;
   double* sat_cut;
   double* press_cut;

   sat_cut = new double[nx_og * ny_og * nz_og];
   sat = new double[nx * ny * nz];
   press = new double[nx * ny * nz];
   press_cut = new double[nx_og * ny_og * nz_og];
   id = new int[nx * ny * nz];
   id_cut = new int[nx_og * ny_og * nz_og]; 

   for (int i = 0 ; i < nprocessors; i++){
   sprintf(nomfic, "mpi_info.%04i", i);   
   Import_MPI_info(nomfic,out);
   rank[i] = out[0];
   locsize[i * 3] = out[1];
   locsize[i * 3 + 1] = out[2];
   locsize[i * 3 + 2] = out[3];	   
   xll[i] = out[7];
   xul[i] = out[8];
   yll[i] = out[9];
   yul[i] = out[10];	   
   zll[i] = out[11];
   zul[i] = out[12];	   
   }
	
   int nindexed = 100;
   int* indexed;
   indexed = new int[nindexed];
	
   for (int i = 0 ; i <= nprocessors; i++){
   Import_id("id_",i,locsize[i * 3],locsize[i * 3 + 1],locsize[i * 3 + 2],id);
   }

   frame_change_id(id,id_cut,nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift,nz_dn_shift,nz_up_shift);

   Export_id("id_ordered_tot",id,(nx * ny * nz));
   Export_id("id_ordered",id_cut,(nx_og * ny_og * nz_og));
   NP = volfrac_calculation(id_cut, nx_og, ny_og, nz_og);	
   NS = (nx_og * ny_og * nz_og) - NP;
   data_mod = new double[NP];

   //save_vtk_ascii("str_tot.vtk", id, nx, ny, nz);
   //save_vtk_ascii("str_cut.vtk", id_cut, nx_og, ny_og, nz_og);

   for (int j = 0 ; j < nindexed ; j++){
   for (int i = 0 ; i < nprocessors; i++){	
   cout<<"reading processor: "<<i<<" for sample: "<<j<<endl;	
   Import_sat("sat_",j,i,locsize[i * 3],locsize[i * 3 + 1],locsize[i * 3 + 2],sat);
   Import_fluid_data("f_",j,i,locsize[i * 3],locsize[i * 3 + 1],locsize[i * 3 + 2],press);
   }

   frame_change(sat,sat_cut,nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift,nz_dn_shift,nz_up_shift);
   Export("rho_",j,sat_cut,(nx_og * ny_og * nz_og));
   // Export("rho_tot_",j,sat,(nx * ny * nz));

   frame_change(press,press_cut,nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift,nz_dn_shift,nz_up_shift);
   Export("p_",j,press_cut,(nx_og * ny_og * nz_og));
   // Export("p_tot_",j,press,(nx * ny * nz));
	   	   
   }

   cout<<"done with analysis"<<endl;

	delete [] data_mod;
	delete [] locsize;
	delete [] out;
	delete [] rank;
	delete [] xll;
	delete [] xul;
	delete [] yll;
	delete [] yul;
	delete [] zll;
	delete [] zul;
	delete [] sat;
	delete [] press;
	delete [] press_cut;
	delete [] id;
	delete [] indexed;
	delete [] id_cut;
	delete [] sat_cut;

   return 0;
}
