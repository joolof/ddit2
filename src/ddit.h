// ----------------------------------------------------------------------------
// Define some constants
// ----------------------------------------------------------------------------
#define HH 6.6262e-27
#define CC2 8.98740441e+20
#define AU 1.496e+13
#define MS 1.99e+33
#define PC 3.08572e+18
#define KK 1.3807e-16
#define mxnang 1000
#define nmxx 1000000
#define CXONE Complex(1.0, 0.0)
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
#define threshold 0.01
// ----------------------------------------------------------------------------
// Define some structures first
// ----------------------------------------------------------------------------
typedef enum { false, true } bool;
// ----------------------------------------------------------------------------
// Structure for the parametes of the disk model
// ----------------------------------------------------------------------------
typedef struct  {
	int nr, ng, nt, nwav, inu_peak;
    double r0, pin, pout, mdisk, dpc, grain, tmin, tmax, geom_corr_fact, smin, smax, density, porosity, fcarbon, fmax;
    char comp[128];
	bool isqfile;
} Param ;
// ----------------------------------------------------------------------------
// Structure for the CLA arguments
// ----------------------------------------------------------------------------
typedef struct  {
	int r0, pin, pout, grain, mdisk;
	bool verbose;
} Cla;
// ----------------------------------------------------------------------------
// Structure for the stellar properties
// ----------------------------------------------------------------------------
typedef struct  {
	double wav, nu, nu3, nuhh, lstar, opt_real, opt_imag, femis ;
} Warray ;
// ----------------------------------------------------------------------------
// Structure for the grain size arrays
// ----------------------------------------------------------------------------
typedef struct  {
    double gsize, qabs, ndens;
} Sarray ;
// ----------------------------------------------------------------------------
// Structure for the radius-temperature array
// ----------------------------------------------------------------------------
typedef struct  {
    double btemp_grid, rt;
} RTarray ;
// ----------------------------------------------------------------------------
void get_input_arguments(int argc, char *argv[], char *star[], Cla *cla, Param *param);
void read_parameters(char *star, Param *param);
int get_numlines(char *filename);
void get_nwav(char *star, Param* param);
void update_parameters(char *argv[], char * star, Param *param, Cla *cla);
void check_validity(Param *param);
void read_stellar(char *star, Warray *warray, Param *param);
void get_optical_cst(Param *param, Warray *warray);
void linint_log_up(double *x, double *y, int n, double *ix, double *iy, int in);
double get_min(double *arr, int len);
double get_max(double *arr, int len);
void bhmie(double x,fcomplex cxref, int nang, fcomplex cxs1[], fcomplex cxs2[], double *qext, double *qsca, float *qback, float *gsca);
void printProgress(double percentage);
void compute_opacity(char *star, Param *param, Sarray *sarray, Warray *warray, RTarray *rtarray);
double get_integrate(double *x, double *fx, int len);
void read_opacity(char *star, Param *param, Sarray *sarray, Warray *warray, RTarray *rtarray);
void get_dust_properties(char *star, Param *param, Sarray *sarray, Warray *warray, RTarray *rtarray, Cla *cla);
void get_sed(char *star, Param *param, Sarray *sarray, Warray *warray, RTarray *rtarray, Cla *cla);
double interpolate_log_down(double *x, double *y, int n, double ix);
int binarySearchDown(double arr[], int h, double x);

