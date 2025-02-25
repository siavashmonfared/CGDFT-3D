#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <set>
#include <string>
#include <sstream>



using namespace std;

double fp, fs;
vector<int>VOL;
// size_t nx, ny, nz;	   
double gridStep; 
int npart;
int Lx,Ly,Lz;
size_t type, nit;
double Tstar, Tcw;
int nxp, nyp, nzp;
double strain, ftol, Rpore; 	
double XC, YC, ZC, obj_tol;
double rho_thresh_ub, rho_thresh_lb;
double rho_thresh;
double wff,wmf,mu,Temp, Boltz,yvar;
double Lx_lt_shift, Lx_rt_shift, Ly_lt_shift, Ly_rt_shift, Lz_up_shift, Lz_dn_shift;

double NP,NS; 
double A11,A12,A13,A14,A15,A16,A17,A18;
double A21,A22,A23,A24,A25,A26,A27,A28;
double A31,A32,A33,A34,A35,A36,A37,A38;
double G11,G12,G13,G21,G22,G23,G31,G32;
double G33,G41,G42,G43,G51,G52,G53,G61;
double G62,G63,G71,G72,G73,G81,G82,G83;
const double pi = 3.1415926535897; 

int ipx(int x, int y, int z, int nnx , int nny){
	return (x + nnx * y + nnx * nny * z);
}




void  Import_data(const string& _name, int num, double tmp[], size_t xl, size_t xu, size_t yl, size_t yu, size_t zl, size_t zu){
	
	double local_var;
	ostringstream str1;
	str1 << num;	
	string app1 = str1.str();
	string result = _name + app1;
	fstream file(result);
	cout<<"just read: "<<result<<endl;
	size_t Nz = zu;
	size_t Ny = yu;
	size_t Nx = xu;
	
	size_t count = 0;
	for (size_t z = 0 ; z < Nz ; z++){					
		for (size_t y = 0 ; y < Ny ; y++) {
			for (size_t x = 0 ; x < Nx ; x++) { 
			file >> local_var;
			if(x >= xl && x < xu && y >= yl && y < yu && z >= zl && z < zu){
			tmp[count] = local_var;
			count++;
			}
			}
		}
	}
}


/*
void  Import_press(const string& _name, int num, double tmp[], size_t xl, size_t xu, size_t yl, size_t yu, size_t zl, size_t zu, double npress){
	
	double local_var;
	string result = _name + to_string(num);
	fstream file(result);
	
	size_t count = 0;
	for (size_t z = 0 ; z < nz ; z++){					
		for (size_t y = 0 ; y < ny ; y++) {
			for (size_t x = 0 ; x < nx ; x++) { 
			file >> local_var;
			if(x >= xl && x < xu && y >= yl && y < yu && z >= zl && z < zu){
			tmp[count] = local_var - npress;
			count++;
			}
			}
		}
	}
}
*/

int  Import_id(const string& _name, int tmpID[], int IDP, size_t xl, size_t xu, size_t yl, size_t yu, size_t zl, size_t zu){

	int tag;
	fstream file(_name);

	size_t Nz = zu;
	size_t Ny = yu;
	size_t Nx = xu;
	
	size_t count = 0;
	int NoPo = 0;
	for (size_t z = 0 ; z < Nz ; z++){					
		for (size_t y = 0 ; y < Ny ; y++) {
			for (size_t x = 0 ; x < Nx ; x++) { 
			file >> tag;
			if(x >= xl && x < xu && y >= yl && y < yu && z >= zl && z < zu){
			tmpID[count] = tag;
			count++;
			if ( tag == IDP ) NoPo++;
			}
			}
		}
	}
return NoPo;
}

/*
size_t volfrac_calculation(int tid[]) 
{
	cout<<"computing volume fractions"<<endl;
	size_t N = nx * ny * nz;
	double tot = N;
	NP = 0.0;
	NS = 0.0;

	int f;
	for (size_t z = 0 ; z < N ; ++z) { 
	f = tid[z];
	if(f==0){NP++;};
	if(f>0){NS++;};
	}

	fp = NP/tot;
	fs = NS/tot;

	cout<<"fp: "<<fp<<" fs: "<<fs<<endl;
	size_t nsi = NS;
	return nsi;
}
*/

double find_min(const string& _name, double tmp[], size_t tot){

	const char * c = _name.c_str();
	
	double small = tmp[0];
	for (size_t i = 0 ; i < tot ; i++){
	if(tmp[i] < small) small = tmp[i];
	}

	FILE * sortie;
	sortie = fopen(c,"a");
	fprintf(sortie,"%g\n",small);
	fclose(sortie);

	return small;

}

double find_max(const string& _name, double tmp[], size_t tot){
	
   	const char * c = _name.c_str();

	double large = tmp[0];
	for (size_t i = 0 ; i < tot ; i++){
	if(tmp[i] > large) large = tmp[i];
	}

	FILE * sortie;
	sortie = fopen(c,"a");
	fprintf(sortie,"%g\n",large);
	fclose(sortie);

	return large;
}

/*
void histogram(const string& _name, int num, double tmp[], size_t tot, double val_min, double val_max, double binwidth){
	
   string result = _name + to_string(num);
   string resultOUT = result;		
   const char * c = resultOUT.c_str();
	
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

	for (size_t i = 0 ; i < tot ; i++){
    int bucket = (int)floor(tmp[i] / binwidth);
    hist[bucket - (int)floor(val_min/binwidth) ] += 1;
	}	

	FILE * sortie;
	sortie = fopen(c,"w+");
	for(int i = 0 ; i < nbins ; i++) fprintf(sortie,"%g %i\n",intv[i],hist[i]);
	fclose(sortie);
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

	for (size_t i = 0 ; i < tot ; i++){
    int bucket = (int)floor(-tmp[i] / binwidth);
    hist[bucket - (int)floor(-val_max/binwidth)] += 1;
	}
	FILE * sortie;
	sortie = fopen(c,"w+");
	for(int i = 0 ; i < nbins ; i++) fprintf(sortie,"%g %i\n",intv[i],hist[nbins - 1 - i]);
	fclose(sortie);
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

	for (size_t i = 0 ; i < tot ; i++){
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
	fclose(sortie);

	delete [] hist_pos;
	delete [] hist_neg;
	delete [] intv;	
	}
		
}
*/


void write_output(const string& _name, int npart, int nbins, double tbins[], double tcounts[], double rho){

   const char * c = _name.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");	

   for (size_t i = 0 ; i < nbins ; i++){
   double interval = 0.5 * (tbins[i] + tbins[i+1]);
   //double freq = tcounts[i]/npart/(4.0*pi*pow(interval,2)*(tbins[i+1]-tbins[i]))/rho;
   double freq = tcounts[i]/(4.0*pi*pow(interval,2)*(tbins[i+1]-tbins[i]))/rho;
   fprintf(sortie,"%.7e %.7e %.7e %.7e %.7e\n",tbins[i],tbins[i+1], interval, freq,rho);
   }
   fclose(sortie);
}



void compute_cumulants(const string& _name,double lamb, double tmp[], int tid[] , int tot, int IDP){
	
	ostringstream str1;
	str1 << lamb;
    string app1 = str1.str();
	
	string result = _name;
	result.append("_");
	string resultIN = result + app1;
	const char * c = resultIN.c_str();

	vector<double> tmp_storage;
	double m1, m2, m4, count, ExprA, ExprB, ExprC, ExprE;
	m1 = m2 = m4 = count = ExprA = ExprB = ExprC = ExprE = 0.0;
	
	for (int i = 0 ; i < tot ; i++){
	if (tid[i] == IDP){
	count++;
	m1 += tmp[i];
	m2 += pow(tmp[i],2);
	m4 += pow(tmp[i],4);
	tmp_storage.push_back(tmp[i]);
	}
	}
	
	m1 = m1/(count);
	m2 = m2/(count);
	m4 = m4/(count);
	
	for (int i = 0 ; i < tmp_storage.size() ; i++){
	ExprE += tmp_storage[i] - m1;
	ExprA += pow((tmp_storage[i] - m1),2);
	ExprB += pow((tmp_storage[i] - m1),3);
	ExprC += pow((tmp_storage[i] - m1),4);
	}

	double dissus = pow((ExprE / count),2);
	double mean = m1;//first cumulant
	double std = sqrt(ExprA/(count - 1));
	double var = pow(std,2);
	double skew = (ExprB/count) / pow(sqrt(ExprA/count),3);
	double kurt = (ExprC/count) / pow(ExprA/count,2);
	
	FILE * sortie;
	sortie = fopen(c,"a");
	fprintf(sortie,"%g %g %g %g %g %g %g\n",mean,var,skew,kurt,m2,m4,dissus);
	fclose(sortie);	
}



int Mask_id(int tid[] , int tid_mask[] , int tot , int IDP , int IDS, int nx, int ny, int nz){
	
	bool xn,yn,zn,xp,yp,zp;
	int* n;
	n = new int[6];
	for (int i = 0 ; i < 6 ; i++) n[i] = 0;
	int count = 0;

 	for (int z = 0 ; z < nz ; z++) {					
	for (int y = 0 ; y < ny ; y++) {
	for (int x = 0 ; x < nx ; x++) {

	xn = (x == nx - 1) ? false : true;
	yn = (y == ny - 1) ? false : true;
	zn = (z == nz - 1) ? false : true;
	xp = (x == 0) ? false : true;
	yp = (y == 0) ? false : true;
	zp = (z == 0) ? false : true;

	int p1 = tid[ipx(x,y,z,nx,ny)];

	if (xn) n[0] = tid[ipx(x+1,y,z,nx,ny)];
	if (yn) n[1] = tid[ipx(x,y+1,z,nx,ny)];
	if (zn) n[2] = tid[ipx(x,y,z+1,nx,ny)];
		
	if (xp) n[3] = tid[ipx(x-1,y,z,nx,ny)];
	if (yp) n[4] = tid[ipx(x,y-1,z,nx,ny)];
	if (zp) n[5] = tid[ipx(x,y,z-1,nx,ny)];

	int sn = 0;
	for (int i = 0 ; i < 6 ; i++) sn += n[i];
	if (p1 == IDP && sn > 0){
	tid_mask[ipx(x,y,z,nx,ny)] = IDS;
	count++;
	}
		
	}
	}
	}
		
	delete [] n;
	return count;
}



void  Export_id(const string& _name, int num, int data[], size_t tot){
	
	
	ostringstream str1;
	str1 << num;
    string app1 = str1.str();
	
	string result = _name;
	result.append("_");
	result = result + app1;
	const char * c = result.c_str();

   FILE * sortie;
   sortie = fopen(c, "w+");	
   for (size_t k = 0 ; k < (tot) ; k++) fprintf(sortie,"%i\n",data[k]);
   fclose(sortie);
}


void  Export_coordinates(const string& _name, int npart, double xc[], double yc[], double zc[], double rad){
	
	
   const char * c = _name.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");	
   for (int k = 0 ; k < (npart) ; k++) fprintf(sortie,"%g %g %g %g\n",xc[k],yc[k],zc[k],rad);
   fclose(sortie);
}


void coarse_graining_func(const string& _name,int num,double lambda,int id[],double data[],double Lx,double Ly,double Lz,int IDP,double gridStep)	{

	int nvx = Lx / lambda + 1;
	int nvy = Ly / lambda + 1;
	int nvz = Lz / lambda + 1;
	int nsites = lambda / gridStep + 1;
	
	int nxn = Lx / gridStep + 1;
	int nyn = Ly / gridStep + 1;
	int nzn = Lz / gridStep + 1;

	int* LengthX;
	int* LengthY;
	int* LengthZ;
	int* intX;
	int* intY;
	int* intZ;
	LengthX = new int[nvx];
	LengthY = new int[nvy];
	LengthZ = new int[nvz];
	intX = new int[nsites];
	intY = new int[nsites];
	intZ = new int[nsites];
	int lbx, ubx, lby, uby, lbz, ubz;
		
	double* val_data;
	val_data = new double[nsites * nsites * nsites];

	int* id_data;
	id_data = new int[nsites * nsites * nsites];		
		
	for (int i = 0 ; i < nvx ; i++) LengthX[i] = i * (nsites - 1);
	for (int i = 0 ; i < nvy ; i++) LengthY[i] = i * (nsites - 1);
	for (int i = 0 ; i < nvz ; i++) LengthZ[i] = i * (nsites - 1);
	
	
	int vol_count = 0;
	for (int i = 0 ; i < nvz-1 ; i++){
	lbz =  LengthZ[i];
	ubz =  LengthZ[i + 1];
	for (int j = 0 ; j < nvy-1 ; j++){
	lby =  LengthY[j];
	uby =  LengthY[j + 1];
	for (int k = 0 ; k < nvx-1 ; k++){
	lbx =  LengthX[k];
	ubx =  LengthX[k + 1];

	for (int ii = 0 ; ii < nsites ; ii++){
	intX[ii] = lbx + ii;	
	intY[ii] = lby + ii;	
	intZ[ii] = lbz + ii;
	}

	size_t count = 0;
	for (int zz = 0 ; zz < nsites ; zz++){
	for (int yy = 0 ; yy < nsites ; yy++){
	for (int xx = 0 ; xx < nsites ; xx++){	
	size_t index = ipx(intX[xx],intY[yy],intZ[zz],nxn,nyn);
	val_data[count] = data[index];
	id_data[count] = id[index];
	count++;
	}
	}
	}
	
	ostringstream str1,str2;
    str1 << num ;
    string app1 = str1.str();
	string result = _name + app1;
	compute_cumulants(result,lambda,val_data,id_data,pow(nsites,3),IDP);
	vol_count++;
	
	}
	}
	}
			
	delete [] LengthX;
	delete [] LengthY;
	delete [] LengthZ;
	delete [] intX;
	delete [] intY;
	delete [] intZ;
	delete [] val_data;
	delete [] id_data;

}



int readSimData(const string& _name, double rh, int TOT){
	
	int a1;
	double a2,a3,a4,a5,a6;
	int index_out;
	ifstream file(_name);

	for (int k = 0 ; k < TOT ; ++k){
	file >> a2 >> a3 >> a4 >> a5 >> a6 >> a1;
	if(a5 == rh) index_out = k;
	}

	return index_out;
}

void  voronoi_volume_frac(const string& _name, int npart, double rad, double fi[] ){
	int ipart;
	double xc, yc, zc, cvol;
	fstream file(_name);
	for (int i = 0 ; i < npart ; i++){
	file >> ipart >> xc >> yc >> zc >> cvol;	
	fi[i] = ( 4.0 * M_PI * pow(rad,3)  / 3.0 ) / cvol;
	}
}


	
void  Import_data_cell( const string& _name, int cpart, int idpart[], double gs, double dist[], double xpar[], double ypar[], double zpar[] ){
	int ipart, index;
	double d, xc, yc, zc;
	fstream file(_name);
	for (int i = 0 ; i < cpart ; i++){
	file >> ipart >> index >> d >> xc >> yc >> zc;	
	dist[i] = d;
	idpart[i] = ipart;
	xpar[ipart - 1] = xc;
	ypar[ipart - 1] = yc;
	zpar[ipart - 1] = zc;
	}
}

void  isolate_data(double tdata[], int tdata_size , double isodata[], int tid[], int id_flag){

	int j = 0;
	for (int i = 0 ; i < tdata_size ; i++){
	if(tid[i] == id_flag){
	isodata[j] = tdata[i];
	j++;
	}	
	}
}


void compute_statistics(const string& _name, vector<double> &vect){
	
	const char * cA = _name.c_str();
	double m1, m2, count, ExprA, ExprB, ExprC;

	m1 = m2 = count = ExprA = ExprB = ExprC = 0.0;
	
	for (int i = 0 ; i < vect.size() ; i++){
	count++;
	m1 += vect[i];
	m2 += pow(vect[i],2);
	}
		
	m1 = m1/(count);
	m2 = m2/(count);
	
	for (int i = 0 ; i < vect.size() ; i++){
	ExprA += pow((vect[i] - m1),2);
	ExprB += pow((vect[i] - m1),3);
	ExprC += pow((vect[i] - m1),4);
	}

	double Mean = m1;//first cumulant
	double Std = sqrt(ExprA/(count - 1));
	//double Var = m2 - pow(m1,2);//second cumulant
	double Var = pow(Std,2);
	double Skew = (ExprB/count) / pow(sqrt(ExprA/count),3);
	double Kurt = (ExprC/count) / pow(ExprA/count,2);
	
	FILE * sortie;
	sortie = fopen(cA,"a");
	if (vect.size() > 0.0)  fprintf(sortie,"%g %g %g %g\n",Mean, Var, Skew, Kurt);
	if (vect.size() == 0.0) fprintf(sortie,"0 0 0 0\n");
	fclose(sortie);
}






void create_xsection(const string& _name, int num, double data[], int tid[], int idp, int ids, int Nz, int Ny, int Nx, int zcut, double flag){

   ostringstream str1;
   str1 << num;
   string app1 = str1.str();
	
   string result = _name + app1;		
   const char * c = result.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");	

   int count = 0;
   for (int z = 0 ; z < Nz ; z++){
	for (int y = 0 ; y < Ny ; y++){ 
	 for (int x = 0 ; x < Nx ; x++){ 
		 if ( z == zcut && tid[count] == idp) fprintf(sortie,"%g\n",data[count]);
		 if ( z == zcut && tid[count] == ids) fprintf(sortie,"%g\n",flag);
	 count++;
	 }
	}
   }
	
	
 fclose(sortie);
	
}



void create_xsection_id(const string& _name, int num, int tid[], int Nz, int Ny, int Nx, int zcut){

   ostringstream str1;
   str1 << num;
   string app1 = str1.str();
   string result = _name + app1;		
   const char * c = result.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");	

   int count = 0;
   for (int z = 0 ; z < Nz ; z++){
   for (int y = 0 ; y < Ny ; y++){ 
   for (int x = 0 ; x < Nx ; x++){ 
   if ( z == zcut ) fprintf(sortie,"%i\n",tid[count]);
   count++;
   }
   }
   }
	
 fclose(sortie);
	
}


void  findGLsites(const string& _name, int tid[], double tdata[], int N, int IDP, double thresh){
	
   int ll = 0;
   int gg = 0;
   int pp = 0;
   for ( int i = 0 ; i < N ; i++){
   if (tid[i] == IDP && tdata[i] >= thresh) ll++;
   if (tid[i] == IDP && tdata[i] < thresh) gg++;
   if (tid[i] == IDP) pp++;
   }
	
   const char * c = _name.c_str();
   FILE * sortie;
   sortie = fopen(c, "a");	
   fprintf(sortie,"%i %i %i\n",pp,ll,gg);
   fclose(sortie);
}

void readSimParam(){
	
    const char * name = "input_12102018.txt";
    ifstream file(name);
    file >> Lx >> Ly >> Lz;
	file >> nxp >> nyp >> nzp;
    file >> XC >> YC >> ZC >> Rpore;
	file >> type >> strain >> nit >> ftol;
	file >> Tstar >> gridStep >> Tcw >> yvar;
	file >> obj_tol >> npart;
	file >> rho_thresh_ub >> rho_thresh_lb;
	file >> Lx_lt_shift >> Lx_rt_shift;
	file >> Ly_lt_shift >> Ly_rt_shift;
	file >> Lz_dn_shift >> Lz_up_shift;
	// file >> rho_thresh;
}


void  Import_sres(const string& _name, double tdata[], int N){
	double msxx,msyy,mszz,msxy,msxz,msyz;
	fstream file(_name);
	for (int i = 0 ; i < N ; i++){					
	file >> msxx >> msyy >> mszz >> msxy >> msxz >> msyz;
	tdata[i] = (1.0/3.0)*(msxx+msyy+mszz);
 	}
}


void resize_data_for_cg(int id[], int idn[], double data[], double datan[], double Lx, double Ly, double Lz, double lambda, double gridStep){
	
	double lxp = lambda - fmod (Lx,lambda) ;
	double lyp = lambda - fmod (Ly,lambda) ;
	double lzp = lambda - fmod (Lz,lambda) ;
		
	int nx = Lx / gridStep + 1;
	int ny = Ly / gridStep + 1;
	int nz = Lz / gridStep + 1;
	
	double Lxn = Lx;
	double Lyn = Ly;
	double Lzn = Lz;

	if (fmod(Lx,lambda) > 0) Lxn = Lxn + lxp;
	if (fmod(Ly,lambda) > 0) Lyn = Lyn + lyp;
	if (fmod(Lz,lambda) > 0) Lzn = Lzn + lzp;
	
	int nxn = Lxn / gridStep + 1;
	int nyn = Lyn / gridStep + 1;
	int nzn = Lzn / gridStep + 1;
			
	int xi,yi,zi,idxRaw,idxNew;
	for (int z = 0 ; z < nzn ; z++){
	for (int y = 0 ; y < nyn ; y++){
	for (int x = 0 ; x < nxn ; x++){
	if (z < nz) zi = z;
	if (y < ny) yi = y;
	if (x < nx) xi = x;
	if (z >= nz) zi = z - nz;
	if (y >= ny) yi = y - ny;
	if (x >= nx) xi = x - nx;
	idxRaw = ipx(xi,yi,zi,nx,ny);
	idxNew = ipx(x,y,z,nxn,nyn);		
	if (z < nz && y < ny) { 
	datan[idxNew] = data[idxRaw];	
	idn[idxNew] = id[idxRaw];	
	}
	else { 
	datan[idxNew] = 0.;
	idn[idxNew] = 1.;
	}
	}
	}
	}

}

void save_vtk_ascii_3D_id(int tid[] , int num, int nnx, int nny, int nnz){
	
    FILE * sortie;
    char nomfic[256];
    sprintf(nomfic, "id%04i.vtk", num);
    sortie = fopen(nomfic, "w");
    fprintf(sortie, "# vtk DataFile Version 2.0\n");
    fprintf(sortie, "Sortie domaine LB+LINK\n");
    fprintf(sortie, "ASCII\n");
    fprintf(sortie, "DATASET RECTILINEAR_GRID\n");
    fprintf(sortie, "DIMENSIONS %i %i %i\n", nnx, nny, nnz);
    
    fprintf(sortie, "X_COORDINATES %i float\n", nnx);
    for (size_t i = 0; i <  nnx; i++) {
        fprintf(sortie, "%.4e ", i * gridStep);
    }
    fprintf(sortie, "\n");
    
    fprintf(sortie, "Y_COORDINATES %i float\n", nny);
    for (size_t i = 0; i < nny; i++) {
        fprintf(sortie, "%.4e ", i * gridStep);
    }
    fprintf(sortie, "\n");
    
    fprintf(sortie, "Z_COORDINATES %i float\n", nnz);
    for (size_t i = 0; i < nnz ; i++) {
        fprintf(sortie, "%.4e ", i * gridStep);
    }
    fprintf(sortie, "\n");
    
    fprintf(sortie, "POINT_DATA %i\n", nnx * nny * nnz);

    size_t ik = 0;
    fprintf(sortie, "SCALARS MatterId int 1\n");
    fprintf(sortie, "LOOKUP_TABLE default\n");
    for (size_t i = 0 ; i < (nnx * nny * nnz) ; i++) fprintf(sortie, "%i\n", tid[i]);
    
    fprintf(sortie, "VECTORS Displacement float\n");
    for (size_t i = 0 ; i < (nnx * nny * nnz) ; i++) fprintf(sortie, "%.4e %.4e %.4e\n", 0.,0.,0.);
    
    fprintf(sortie, "VECTORS StressDiag float\n");
    for (size_t i = 0 ; i < (nnx * nny * nnz) ; i++) fprintf(sortie, "%.4e %.4e %.4e\n", 0.,0.,0.);
    
    fclose(sortie);
}

// Writes the 3D array 'data' to a VTK ASCII file "fname".
// The data is assumed to be shaped (nz, nx, ny) in row-major layout:
//    data[z*(nx*ny) + x*ny + y]
// matching how the Python version loops z, x, y.
//
// For compatibility with the Python snippet:
//   - "DIMENSIONS" is written as (ny, nx, nz).
//   - The data is written in nested loops over z, x, y.
//
// Adjust as needed if your indexing or dimension order is different.

void save_vtk_ascii_3D(double data[] , int num, int nx, int ny, int nz){
	
    FILE * f;
    char nomfic[256];
    sprintf(nomfic, "vol%04i.vtk", num);
    f = fopen(nomfic, "w");

    // VTK header
    std::fprintf(f, "# vtk DataFile Version 2.0\n");
    std::fprintf(f, "volume example\n");
    std::fprintf(f, "ASCII\n");
    std::fprintf(f, "DATASET STRUCTURED_POINTS\n");

    // Note: Python snippet did "DIMENSIONS %d %d %d" with (Y,X,Z).
    // We'll replicate that exact ordering.
    std::fprintf(f, "DIMENSIONS %d %d %d\n", nx, ny, nz);
    std::fprintf(f, "ASPECT_RATIO 1 1 1\n");
    std::fprintf(f, "ORIGIN 0 0 0\n");

    // Write scalar data
    int totalPoints = nx * ny * nz;
    std::fprintf(f, "POINT_DATA %d\n", totalPoints);
    std::fprintf(f, "SCALARS volume_scalars double 1\n");
    std::fprintf(f, "LOOKUP_TABLE default\n");

    // Loop in z, x, y order to match how the Python code does:
    //   for z in [0..nz):
    //       for x in [0..nx):
    //           for y in [0..ny):
    //               print(data[z, x, y])
    /*
    for (int z = 0; z < nz; z++) {
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                double val = data[z * (nx * ny) + x * ny + y];
                std::fprintf(f, "%.6f\n", val);
            }
        }
    }
    */
    for (size_t i = 0 ; i < (nx * ny * nz) ; i++) fprintf(f, "%.4f\n", data[i]);

    std::fclose(f);
}

/*
void save_vtk_ascii_3D(double tdata[] , int num, int nnx, int nny, int nnz){
	
    FILE * sortie;
    char nomfic[256];
    sprintf(nomfic, "vol%04i.vtk", num);
    sortie = fopen(nomfic, "w");
    fprintf(sortie, "# vtk DataFile Version 2.0\n");
    fprintf(sortie, "Sortie domaine LB+LINK\n");
    fprintf(sortie, "ASCII\n");
    fprintf(sortie, "DATASET RECTILINEAR_GRID\n");
    fprintf(sortie, "DIMENSIONS %i %i %i\n", nnx, nny, nnz);
    
    fprintf(sortie, "X_COORDINATES %i float\n", nnx);
    for (size_t i = 0; i <  nnx; i++) {
        fprintf(sortie, "%.4e ", i * gridStep);
    }
    fprintf(sortie, "\n");
    
    fprintf(sortie, "Y_COORDINATES %i float\n", nny);
    for (size_t i = 0; i < nny; i++) {
        fprintf(sortie, "%.4e ", i * gridStep);
    }
    fprintf(sortie, "\n");
    
    fprintf(sortie, "Z_COORDINATES %i float\n", nnz);
    for (size_t i = 0; i < nnz ; i++) {
        fprintf(sortie, "%.4e ", i * gridStep);
    }
    fprintf(sortie, "\n");
    
    fprintf(sortie, "POINT_DATA %i\n", nnx * nny * nnz);

    size_t ik = 0;
    fprintf(sortie, "SCALARS MatterId int 1\n");
    fprintf(sortie, "LOOKUP_TABLE default\n");
    for (size_t i = 0 ; i < (nnx * nny * nnz) ; i++) fprintf(sortie, "%i\n", 1);
    
    fprintf(sortie, "VECTORS Displacement float\n");
    for (size_t i = 0 ; i < (nnx * nny * nnz) ; i++) fprintf(sortie, "%.4e %.4e %.4e\n", tdata[i],0.,0.);
    
    fprintf(sortie, "VECTORS StressDiag float\n");
    for (size_t i = 0 ; i < (nnx * nny * nnz) ; i++) fprintf(sortie, "%.4e %.4e %.4e\n", 0.,0.,0.);
    
    fclose(sortie);
}
*/



int main (){
	
	readSimParam();
	

	int IDP = 0;
	int IDS = 1;
	int nsim = 100;

	// #define NoLamb 7
	// const double lamb[NoLamb] = {11,12,13,14,15,16,17};//This is length not number of nodes
	
	#define NoLamb 8
	const double lamb[NoLamb] = {1,2,3,4,5,6,7,8};//This is length not number of nodes

	int nx = Lx / gridStep + 1;
	int ny = Ly / gridStep + 1;
	int nz = Lz / gridStep + 1;
		
	int* id;
	double* sat;

   	id 	=  new int[nx * ny * nz];
   	sat =  new double[nx * ny * nz];

   	int NoPo_r = Import_id("id_ordered", id, IDP,0,nx,0,ny,0,nz);
	
	// for (int i = 0 ; i < (nx*ny*nz) ; i++) id_mask[i] = id[i];
	// int NoPo_sub = Mask_id(id,id_mask,(nx*ny*nz),IDP,IDS,nx,ny,nz);
	// int NoPo_p = NoPo_r - NoPo_sub;
	// create_xsection_id("id_xsec_", 1, id, nz, ny, nx, 40);
	// create_xsection_id("id_mask_xsec_", 1, id_mask, nz, ny, nx, 40);
	
	for (int i = 0 ; i < nsim ; i++){
	cout<<"importing data :"<<i<<endl;
	// int i = 88;
	Import_data("rho_",i, sat,0,nx,0,ny,0,nz);
	// save_vtk_ascii_3D_id(id,1,nx,ny,nz);	
	save_vtk_ascii_3D(sat,i,nx,ny,nz);	
	// save_vtk_ascii_3D_id(id_cropped,2,nx,(nyf-nyi),(nzf-nzi));	
	// save_vtk_ascii_3D(sat_cropped,2,nx,(nyf-nyi),(nzf-nzi));	
	
	for (size_t jj = 0 ; jj < NoLamb ; jj++){
    
	cout<<"reading sample: "<<i<<" cg length scale: "<<lamb[jj]<<endl;		
	double lambda = lamb[jj];
		
	double lxp = lambda - fmod (Lx,lambda) ;
	double lyp = lambda - fmod (Ly,lambda) ;
	double lzp = lambda - fmod (Lz,lambda) ;
	
	double Lxn = Lx;
	double Lyn = Ly;
	double Lzn = Lz;

	if (fmod(Lx,lambda) > 0) Lxn = Lxn + lxp;
	if (fmod(Ly,lambda) > 0) Lyn = Lyn + lyp;
	if (fmod(Lz,lambda) > 0) Lzn = Lzn + lzp;
		
	int nxn = Lxn / gridStep + 1;
	int nyn = Lyn / gridStep + 1;
	int nzn = Lzn / gridStep + 1;
				
	double* datan;
	int* idn;	
		
	datan = 		new double[nxn * nyn * nzn];
	idn = 		new int[nxn * nyn * nzn];		
		
	resize_data_for_cg(id,idn,sat,datan,Lx,Ly,Lz,lambda,gridStep);
	coarse_graining_func("r_cg_rec_",i,lambda,idn,datan,Lxn,Lyn,Lzn,IDP,gridStep);
	// save_vtk_ascii_3D_id(idn,3,nxn,nyn,nzn);	
	// save_vtk_ascii_3D(datan,3,nxn,nyn,nzn);
		
	delete [] datan;
	delete [] idn;
   	}
	}

	delete [] id;
	delete [] sat;
	
	cout<<"computation done."<<endl;
	return 0;
}


	
