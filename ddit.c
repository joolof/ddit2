#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "complex.h"
#include "ddit.h"
#include "nrutil.h"
/* 
---------------------------------------------------------
Main program, call all the functions in order
---------------------------------------------------------
gcc -o ddit ddit.c utils.c complex.c nrutil.c -lm
*/
int main(int argc, char *argv[])
{
	char * star ;
	double seconds ;
	//  ------------------------------------
    clock_t start = clock();
	//  ------------------------------------
	// Initialise all the structures
	//  ------------------------------------
	Param param = {0, 0, 0, 0, 0, 0, 0., 0., 0., 0., 0., 0., 0., 3., 2000., 0., 0.1, 1000., 3.5, "", false};
	Cla cla = {0, 0, 0, 0, 0, 0, false};
    Sarray * sarray ; // structure for grain sizes dependant arrays
    RTarray * rtarray ; // structure for radius-temperature arrays
	Warray * warray ; // structure for the wavelength dependent array
	//  ------------------------------------
	//  Get the name of the parameter file
	//  ------------------------------------
	get_input_arguments(argc, argv, &star, &cla, &param);
	//  ------------------------------------
	//  Read the parameter file
	//  ------------------------------------
	read_parameters(star, &param) ;
	get_nwav(star, &param);
	//  ------------------------------------
	//  Change some paraneters if they were passed as CLA
	//  ------------------------------------
	//printf("%f\n", param.smin);
	update_parameters(argv, star, &param, &cla);
	//printf("%f\n", param.smin);
	check_validity(&param);
	//  ------------------------------------
	//  Initialize the structure and read the stellar properties
	//  ------------------------------------
    sarray = malloc ((param.nwav * param.ng) * sizeof(Sarray));
    rtarray = malloc ((param.nt * param.ng) * sizeof(RTarray));
	warray = malloc (param.nwav * sizeof(Warray));
	read_stellar(star, warray, &param);
    // ------------------------------------
	// Get the dust properties
    // ------------------------------------
    get_dust_properties(star, &param, sarray, warray, rtarray, &cla);
    // ------------------------------------
	// Do the thing
    // ------------------------------------
	get_sed(star, &param, sarray, warray, rtarray);
    // ------------------------------------
    // Finished
    // ------------------------------------
    clock_t end = clock();
    seconds = (float)(end - start) / CLOCKS_PER_SEC;
    if (cla.verbose)
    {
        printf("--------------------------------------------------------------------------------\n");
        printf("CPU time\t\t= %.2f\t\t[s]\n",seconds);
        printf("--------------------------------------------------------------------------------\n\n");
    }
    free(warray);
    free(sarray);
    free(rtarray);
    return 0;
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Ok, let's do the thing
---------------------------------------------------------
---------------------------------------------------------
*/
void get_sed(char *star, Param *param, Sarray *sarray, Warray *warray, RTarray *rtarray)
{
	int ir, ig, iwav, iz, it;
	char junk[1000];
	double ri[param->nr+1], zi[param->nz+1], rtemp[param->nt], btemp[param->nt];
	double dist = 0., rdens = 0., rin = 0., rout = 0., rmid = 0., height = 0., vol = 0.;
	double zmid = 0., tgrain = 0., tmp = 0., thermc = 0., mdust = 0., Tmax = -1.;
	// -------------------------------------------	
    for (it=0 ; it < param->nt ; it++) btemp[it] = rtarray[it].btemp_grid ;
	for (iwav = 0; iwav < param->nwav; iwav++) warray[iwav].femis = 0.;
	rin = pow(threshold, 1.e0/param->pin) * param->r0 ;
	rout = pow(threshold, 1.e0/param->pout) * param->r0 ;
	// -------------------------------------------	
	for (ir=0 ;ir < param->nr + 1; ir++) 
	{
		ri[ir] = rin * pow((rout / rin),((ir * 1.0) / (param->nr * 1.0)));
	}
	// -------------------------------------------	
	for (ir=0 ;ir < param->nr; ir++) 
	{
		rmid = sqrt(ri[ir]*ri[ir+1]);
		for (iz = 0; iz < param->nz + 1; iz++)
		{
			height = rmid * param->opang;
			zi[iz] = -5.*height + 10.*height*((iz*1.)/(param->nz*1.));
		}
		for (iz = 0; iz < param->nz; iz++)
		{
			zmid = (zi[iz] + zi[iz+1])/2.e0;
			dist = sqrt(zmid * zmid + rmid * rmid);
			rdens = exp(-1.e0 * zmid * zmid / (2. * height * height)) / sqrt(pow(rmid / param->r0, -2.e0 * param->pin) + pow(rmid/param->r0, -2.e0 * param->pout));
			vol = M_PI * ((ri[ir+1] * ri[ir+1]) - (ri[ir] * ri[ir])) * (zi[iz+1] - zi[iz]);
			for (ig = 0; ig < param->ng; ig++)
			{
				for (it = 0; it < param->nt ; it++)
				{
					rtemp[it] = rtarray[it * param->ng + ig].rt;
				}
				tgrain = interpolate_log_down(rtemp, btemp, param->nt, dist);
				if (tgrain > Tmax) Tmax = tgrain;
				for (iwav = 0; iwav < param->nwav; iwav++)
				{
					tmp = exp(warray[iwav].nuhh / tgrain);
					thermc = 2.0 * HH * M_PI * warray[iwav].nu3 / CC2 / (tmp - 1.) ;
					thermc *= 4.0 * M_PI *  sarray[ig].gsize *sarray[ig].gsize * sarray[iwav * param->ng + ig].qabs;
					tmp = vol * rdens * sarray[ig].ndens;
					warray[iwav].femis += (tmp * thermc);
				}
				mdust += (tmp * 4.0/3.0 * M_PI * param->density * sarray[ig].gsize * sarray[ig].gsize * sarray[ig].gsize);
			}
		}

	}
	printf("%f\n", Tmax);
  	sprintf(junk, "%s/SED.dat",star);
	FILE *fichier = fopen(junk, "w" );
	fprintf(fichier,"Wave Lstar Disk\n");
	for (iwav = 0; iwav < param->nwav; iwav++)
	{
		warray[iwav].femis = warray[iwav].femis / mdust * param->mdisk * param->geom_corr_fact;
		fprintf(fichier,"%f %8.10e %8.10e\n", warray[iwav].wav, warray[iwav].lstar, warray[iwav].femis);
	}
	fclose(fichier);
	

}
/* 
---------------------------------------------------------
---------------------------------------------------------
Read the parameter file
---------------------------------------------------------
---------------------------------------------------------
*/
void read_parameters(char *star, Param *param)
{
	char junk[1000] ;
	sprintf(junk, "%s/%s.ddit",star,star);
	FILE *fparam = fopen(junk, "r" );
	if (fparam == 0)
	{
		printf("Could not open the param file (%s/%s.ddit).\n", star, star);
		printf("--------------------------------------------------------------------------------\n");
		exit(0);
	}
    else 
    {
		while (fgets(junk, 1000, fparam) != 0)
		{
			sscanf(junk, "dpc %lf", &param->dpc) ;
			sscanf(junk, "r0 %lf", &param->r0) ;
			sscanf(junk, "pin %lf", &param->pin) ;
			sscanf(junk, "pout %lf", &param->pout) ;
			sscanf(junk, "grain %lf", &param->grain) ;
			sscanf(junk, "smin %lf", &param->smin) ;
			sscanf(junk, "smax %lf", &param->smax) ;
			sscanf(junk, "mdisk %lf", &param->mdisk) ;
			sscanf(junk, "opang %lf", &param->opang) ;
			sscanf(junk, "nr %d", &param->nr) ;
			sscanf(junk, "ng %d", &param->ng) ;
			sscanf(junk, "nz %d", &param->nz) ;
			sscanf(junk, "nt %d", &param->nt) ;
			sscanf(junk, "comp %s", param->comp) ;
		}
		fclose(fparam);
    }
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Check the number of input arguments, and get the 
name of the parameter file. 
---------------------------------------------------------
---------------------------------------------------------
*/
void get_input_arguments(int argc, char *argv[], char *star[], Cla *cla, Param *param)
{
	int i;
	if ((argc < 2) || ((argc ==2) && ((strcmp(argv[1],"-help") == 0) || (strcmp(argv[1],"--help") == 0) || (strcmp(argv[1],"--version") == 0) )))
	{
		if ((argc < 2) || ((argc ==2) && ((strcmp(argv[1],"-help") == 0) || (strcmp(argv[1],"--help") == 0))))
		{
			printf("--------------------------------------------------------------------------------\n");
			printf("For the program to run, you need:\n");
			printf("    + a directory called \"starname\"\n");
			printf("    + a parameter file in this directory \"starname.ddit\"\n");
			printf("                         (see example & documentation)\n");
			printf("    + the directory should additionally contain:\n");
			printf("        - a stellar model, with wavelength and stellar SED\n");
			printf("        - a dust model, with Qabs and Qsca\n");
			printf(" \n");
			printf("--------------------------------------------------------------------------------\n");
			printf(" \n");
			printf("The easiest way to run the program is:\n");
			printf(" \n");
			printf("> %s starname\n",argv[0]);
			printf(" \n");
			printf("--------------------------------------------------------------------------------\n");
			printf(" \n");
			printf("It is also possible to give a few arguments for the disk model:\n");
			printf(" \n");
			printf("> %s starname -r0 40. -pin 10. -pout -3.0 -grain -3.5 -mdisk 1.e-8\n", argv[0]);
			printf(" \n");
			printf("The possible arguments are:\n");
			printf("    -r0 <value>: reference radius (requires pin and pout to be defined)\n");
			printf("    -pin <value>: power-law index for the inner regions\n");
			printf("    -pout <value>: power-law index for the outer regions\n");
			printf("    -grain <value>: power-law index for the grain size dsitribution\n");
			printf("    -mdisk <value>: dust mass\n");
			printf("    -opang <value>: opening angle of the disk\n");
			printf(" \n");
			printf("To print some information while running the code,\nyou can add the keyword \"/verbose\"\n");
			printf(" \n");
			printf("\n");
			printf("--------------------------------------------------------------------------------\n");
			printf(" \n");
		}
		else
		{
			printf("%s 1.0\n", argv[0]);
		}
		exit(0);
	  }
	  else
	  {
		  * star = argv[1];
	  }
	  if (argc > 2)
	  {
			for (i = 2 ; i < argc ; i++)
			{
				if (strcmp(argv[i],"-r0") == 0) cla->r0 = i+1 ;
				if (strcmp(argv[i],"-pin") == 0) cla->pin = i+1;
				if (strcmp(argv[i],"-pout") == 0) cla->pout = i+1;
				if (strcmp(argv[i],"-grain") == 0) cla->grain = i+1;
				if (strcmp(argv[i],"-mdisk") == 0) cla->mdisk = i+1;
				if (strcmp(argv[i],"/verbose") == 0) cla->verbose = true;
				if (strcmp(argv[i],"-opang") == 0) cla->opang = i+1;
			}
	  }

	  if (cla->verbose)
	  {
			printf(" \n");
			printf("Reading parameter file \"%s/%s.ddit\"\n", * star, * star);
			printf("--------------------------------------------------------------------------------\n");
	  }
}
/*
---------------------------------------------------------
---------------------------------------------------------
Get size of the stellar spectrum 
---------------------------------------------------------
---------------------------------------------------------
*/
void get_nwav(char *star, Param* param)
{
  char filename[1000] ;
    sprintf(filename, "%s/starspec.dat",star) ;
	if( access(filename, F_OK ) == -1 ) 
	{
		printf("The stellar spectrum file (%s) does not seem to exist.\n", filename);
		printf("--------------------------------------------------------------------------------\n");
		exit(0);
	}
    param->nwav = get_numlines(filename) ;
}
int get_numlines(char *filename)
{
	int nlines = 0;
	FILE* myfile = fopen(filename, "r");
	int ch;
	do 
	{
		ch = fgetc(myfile);
		if(ch == '\n') nlines = nlines + 1 ;
	} while (ch != EOF);
	fclose(myfile);
	return nlines;
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Update disk parameters if necessary (from CLA)
---------------------------------------------------------
---------------------------------------------------------
*/
void update_parameters(char *argv[], char * star, Param *param, Cla *cla)
{
    char filename[1000];
	if (cla->r0 > 0) param->r0  = atof(argv[cla->r0]) ;
	if (cla->pin > 0) param->pin  = atof(argv[cla->pin]);
	if (cla->pout > 0) param->pout  = atof(argv[cla->pout]);
	if (cla->grain > 0) param->grain  = atof(argv[cla->grain]) ;
	if (cla->opang > 0) param->opang  = atof(argv[cla->opang]) ;
	if (cla->mdisk > 0) param->mdisk  = atof(argv[cla->mdisk]) ;
	if (param->nz % 2 == 0) param->nz+=1;
    //---------------------------------------------------------
    // Check if the opacities were already calculated
    //---------------------------------------------------------
    sprintf(filename, "%s/%s_%.2f_%.2f_%d.dat",star, param->comp, param->smin, param->smax, param->ng);
    if( access(filename, F_OK ) == -1 ) 
    {
        // it does not exist
        param->isqfile = false ;
    }
    else
    {
        // it does exist
        param->isqfile = true ;
    }
	// -------------------------------------------------------
	// Update some parameters
	// -------------------------------------------------------
	param->smin *= 1.e-4;
	param->smax *= 1.e-4;
	param->r0 *= AU ;
	param->mdisk *= MS ;
	param->dpc *= PC ;
	param->geom_corr_fact = 1.e0 / (4.e0 * M_PI * param->dpc * param->dpc);
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Read the frequency file
---------------------------------------------------------
---------------------------------------------------------
*/
void read_stellar(char *star, Warray *warray, Param *param)
{
	int iwav;
	double nulnu, dum;
	char file[1000];
	sprintf(file, "%s/starspec.dat",star);
	// There is already a check on the file beforehand. So I know the file exits.
	FILE *fichier = fopen(file, "r" );
	for (iwav=0; iwav < param->nwav; iwav++)
	{   
		fscanf(fichier,"%lf %lf", &warray[iwav].wav, &warray[iwav].lstar);
		warray[iwav].nu = 2.9979e14/warray[iwav].wav;
		warray[iwav].nu3 = warray[iwav].nu * warray[iwav].nu * warray[iwav].nu ;
		warray[iwav].nuhh = HH * warray[iwav].nu / KK;
	}
	fclose(fichier);
	param->inu_peak = -1;
	dum = -1.;
	for (iwav=0; iwav < param->nwav; iwav++)
	{   
		nulnu =  warray[iwav].lstar * param->dpc * param->dpc * 4.e0 * M_PI * warray[iwav].nu ;
		if (nulnu > dum)
		{
			param->inu_peak = iwav ;
			dum = nulnu ;
		}
	}
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Check the validity of a few parameters before running the
rest of the program.
---------------------------------------------------------
---------------------------------------------------------
*/
void check_validity(Param *param)
{
	// -------------------------------------------------------
	// Make a few checks on the input parameters.
	// -------------------------------------------------------
	if (param->pin <= 0.)
	{
		printf("\'Pin\' cannot be smaller than 0.\n");
		printf("Quitting.\n");
		printf("--------------------------------------------------------------------------------\n");
		exit(0);						
	}
	if (param->pout >= 0.)
	{
		printf("\'Pout\' cannot be larger than 0.\n");
		printf("Quitting.\n");
		printf("--------------------------------------------------------------------------------\n");
		exit(0);						
	}
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Function to get the optical constants
---------------------------------------------------------
---------------------------------------------------------
*/
void get_optical_cst(Param *param, Warray *warray)
{
    char file[1000];
    register unsigned int iwav, nwav_lnk;
    double wav[param->nwav], real[param->nwav], imag[param->nwav];
    //---------------------------------------------------------
    // Check if the environment variable is set or not
    //---------------------------------------------------------
    if (getenv("DDIT_PATH") == NULL)
    {
        printf("Cannot find the file containing the optical constants.\n");
        printf("Please set the environment variable \"DDIT_PATH\" where the .lnk files are stored.\n");
        printf("--------------------------------------------------------------------------------\n");
        exit(0);
    }
    sprintf(file, "%s/%s.lnk", getenv("DDIT_PATH"), param->comp);
    nwav_lnk = get_numlines(file);
    double wav_lnk[nwav_lnk], opt_rlnk[nwav_lnk], opt_ilnk[nwav_lnk];

    FILE *fichier = fopen(file, "r" );
    if (fichier == NULL)
    {
        printf("Problem to read the optical constants file\n");
        printf("\"%s\"\n",file);
        printf("Quitting\n");
        printf("--------------------------------------------------------------------------------\n");
        exit(0);
    }
    else
    {
        for (iwav=0; iwav < nwav_lnk; iwav++)
        {
            fscanf(fichier,"%lf %lf %lf", &wav_lnk[iwav], &opt_rlnk[iwav], &opt_ilnk[iwav]);
        }
        //---------------------------------------------------------
        // Interpolate at the wavelengths of the photospheric model
        // First, need to check if the min and max wavelengths agree
        //---------------------------------------------------------
        for (iwav=0; iwav < param->nwav; iwav++)
        {
            wav[iwav] = warray[iwav].wav ;
        }
        if (get_min(wav_lnk, nwav_lnk) > get_min(wav, param->nwav))
        {
            printf("The minimum wavelength of the optical constant file should be smaller\n");
            exit(0);
        }
        if (get_max(wav_lnk, nwav_lnk) < get_max(wav, param->nwav))
        {
            printf("The maximum wavelength of the optical constant file should be larger\n");
            exit(0);
        }
        linint_log_up(wav_lnk, opt_rlnk, nwav_lnk, wav, real, param->nwav);
        linint_log_up(wav_lnk, opt_ilnk, nwav_lnk, wav, imag, param->nwav);
        for (iwav = 0; iwav < param->nwav; iwav++)
        {
            warray[iwav].opt_real = real[iwav];
            warray[iwav].opt_imag = imag[iwav];
        }
    }
    fclose(fichier);
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Interpolation function, in log-space, for x-arrays in increasing order
---------------------------------------------------------
---------------------------------------------------------
*/
void linint_log_up(double *x, double *y, int n, double *ix, double *iy, int in)
{
    int i = 0, j = 0;
    double logx[2], logy[2];
    double lix[in];
    double lx[n], ly[n];

    for (i =0 ; i < n ; i++)
    {
        lx[i] = log(x[i]);
        ly[i] = log(y[i]);
    }
    for (i =0 ; i < in ; i++)
    {
        lix[i] = log(ix[i]);
    }

    for (i =0 ; i < in ; i++)
    {
        j = 0;
        while (x[j] < ix[i])
        {
            j++;
        }
    
        if (j == 0)
        {
            logx[0] = lx[j];
            logx[1] = lx[j+1];
            logy[0] = ly[j];
            logy[1] = ly[j+1];
        }
        else
        {
            logx[0] = lx[j-1];
            logx[1] = lx[j];
            logy[0] = ly[j-1];
            logy[1] = ly[j];
        }
        iy[i] = logy[0] + (logy[1] - logy[0]) / ( logx[1] - logx[0]) * (lix[i] - logx[0]);
        iy[i] = exp(iy[i]);
    }
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Get the minimum value of an array of size len
---------------------------------------------------------
---------------------------------------------------------
*/
double get_min(double *arr, int len)
{
    int i;
    double min = arr[0];
    for (i = 0; i < len; i++)
    {   
        if (arr[i] < min)
        {
            min = arr[i] ;
        }
    }
    // printf("%f\n", min);
    return min;
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Get the maximum value of an array of size len
---------------------------------------------------------
---------------------------------------------------------
*/
double get_max(double *arr, int len)
{
    int i;
    double max = arr[0];
    for (i = 0; i < len; i++)
    {   
        if (arr[i] > max)
        {
            max = arr[i] ;
        }
    }
    return max;
}
/* 
---------------------------------------------------------
---------------------------------------------------------
bhmie to compute the opacity
---------------------------------------------------------
---------------------------------------------------------
*/
void  bhmie(double x, fcomplex cxref, int nang, fcomplex cxs1[], fcomplex cxs2[], double *qext, double *qsca, float *qback, float *gsca)
{
/*
     Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
     to calculate scattering and absorption by a homogenous isotropic
     sphere.
  Given:
     X = 2*pi*a/lambda
     REFREL = (complex refr. index of sphere)/(real index of medium)
     NANG = number of angles between 0 and 90 degrees
            (will calculate 2*NANG-1 directions from 0 to 180 deg.)
  Returns:
     S1(1 .. 2*NANG-1) =  (incid. E perp. to scatt. plane,
                                 scatt. E perp. to scatt. plane)
     S2(1 .. 2*NANG-1) =  (incid. E parr. to scatt. plane,
                                 scatt. E parr. to scatt. plane)
     QEXT = C_ext/pi*a**2 = efficiency factor for extinction
     QSCA = C_sca/pi*a**2 = efficiency factor for scattering
     QBACK = 4*pi*(dC_sca/domega)/pi*a**2
           = backscattering efficiency
     GSCA = <cos(theta)> for scattering
 
  Original program taken from Bohren and Huffman (1983), Appendix A
  Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
  in order to compute <cos(theta)>
  This code was translatted to C by P. J. Flatau Feb 1998. The C
  version uses "Numerical Recipes" public domain code for complex
  arithmetics "complex.c" and "nrutil.c" (http://www.nr.com).
*/

/* .. Array Arguments .. */
/*      COMPLEX :: cxs1(2*mxnang-1), cxs2(2*mxnang-1)*/
/* .. Local Scalars ..*/
    fcomplex cxan, cxan1, cxbn, cxbn1, cxxi, cxy, cxxi1;
    // fcomplex cxan, cxan1, cxbn, cxbn1, cxxi, cxxi0, cxy, cxxi1;
    fcomplex cxtemp;
    long double  apsi, apsi1, chi, chi0, chi1, dang, fn, p, pii, rn, t, theta, xstop, ymod;
    // long double  apsi, apsi0, apsi1, chi, chi0, chi1, dang, fn, p, pii, rn, t, theta, xstop, ymod;
    long double  dn, dx, psi, psi0, psi1;
    unsigned int  j, jj, n, nmx, nn, nstop;
    /* .. Local Arrays ..*/
    /*fcomplex cxd[nmxx];*/
    fcomplex (*cxd) = malloc(nmxx * sizeof(fcomplex));

    long double  amu[mxnang], pi[mxnang], pi0[mxnang], pi1[mxnang], tau[mxnang];

    if (nang>mxnang)
    {
        printf(" STOP '***Error: NANG > MXNANG in bhmie");
        free(cxd);
        return;
    }
    pii = 4.E0*atan(1.E0);
    dx = x;
    cxy = Cmul(Complex(x,0.0),cxref);

    /* Series expansion terminated after NSTOP terms */
    xstop = x + 4.E0*pow(x,0.3333) + 2.0;
    nstop = xstop;
    ymod = Cabs(cxy);
    nmx = FMAX(xstop,ymod) + 15;

    if (nmx>nmxx) 
    {
        printf(" x, nmx, nmxx, cxref %lf %i %i  \n ", x, nmx, nmxx);
        printf(" xstop nstop ymod %LF %i %LF \n", xstop, nstop, ymod); 
        printf(" Error: NMX > NMXX= %i \n", nmxx);
        free(cxd);
        return;
    }
      

    dang = .5E0*pii/ (float)(nang-1);
    for (j = 1; j<=nang; j++) 
    {
        theta = (float)(j-1)*dang;
        amu[j] = cos(theta);
    }

    /* Logarithmic derivative D(J) calculated by downward recurrence
    beginning with initial value (0.,0.) at J=NMX */
    cxd[nmx] = Complex(0.E0,0.E0);
    nn = nmx - 1;

    for (n = 1; n<= nn; n++) 
    {
        rn = nmx - n + 1;
        // cxd(nmx-n) = (rn/cxy) - (1.E0/(cxd(nmx-n+1)+rn/cxy))
        cxtemp=Cadd(cxd[nmx-n+1],Cdiv(Complex(rn,0.0),cxy));
        cxtemp=Cdiv(CXONE,cxtemp);
        cxd[nmx-n]=Csub(Cdiv(Complex(rn,0.0),cxy),cxtemp);
    }

    for ( j = 1; j <= nang; j++) 
    {
        pi0[j] = 0.E0;
        pi1[j] = 1.E0;
    }
    nn = 2*nang - 1;
    for(j = 1; j<= nn; j++) 
    {
        cxs1[j] = Complex(0.E0,0.E0);
        cxs2[j] = Complex(0.E0,0.E0);
    }
    /* Riccati-Bessel functions with real argument X
    calculated by upward recurrence */
    psi0 = cos(dx);
    psi1 = sin(dx);
    chi0 = -sin(x);
    chi1 = cos(x);
    // apsi0 = psi0;
    apsi1 = psi1;
    // cxxi0 = Complex(apsi0,-chi0);
    cxxi1 = Complex(apsi1,-chi1);
    *qsca = 0.E0;
    *gsca = 0.E0;

    for ( n = 1; n <= nstop; n++) 
    {
        dn = n;
        rn = n;
        fn = (2.E0*rn+1.E0)/(rn*(rn+1.E0));
        psi = (2.E0*dn-1.E0)*psi1/dx - psi0;
        apsi = psi;
        chi = (2.E0*rn-1.E0)*chi1/x - chi0;
        cxxi = Complex(apsi,-chi);
        /* Store previous values of AN and BN for use
        in computation of g=<cos(theta)> */
        if (n>1) 
        {
          cxan1 = cxan;
          cxbn1 = cxbn;
        }
        /* Compute AN and BN:*/
        // cxan = (cxd(n)/cxref+rn/x)*apsi - apsi1;
        cxan=Cdiv(cxd[n],cxref);
        cxan=Cadd(cxan,Complex(rn/x,0.0));
        cxan=Cmul(cxan,Complex(apsi,0.0));
        cxan=Csub(cxan,Complex(apsi1,0.0));
        // cxan = cxan/((cxd(n)/cxref+rn/x)*cxxi-cxxi1);
        cxtemp=Cdiv(cxd[n],cxref);
        cxtemp=Cadd(cxtemp,Complex(rn/x,0.0));
        cxtemp=Cmul(cxtemp,cxxi);
        cxtemp=Csub(cxtemp,cxxi1);
        cxan=Cdiv(cxan,cxtemp);
        // cxbn = (cxref*cxd(n)+rn/x)*apsi - apsi1;
        cxbn=Cmul(cxref,cxd[n]);
        cxbn=Cadd(cxbn,Complex(rn/x,0.0));
        cxbn=Cmul(cxbn,Complex(apsi,0.0));
        cxbn=Csub(cxbn,Complex(apsi1,0.0));
        // cxbn = cxbn/((cxref*cxd(n)+rn/x)*cxxi-cxxi1);
        cxtemp=Cmul(cxref,cxd[n]);
        cxtemp=Cadd(cxtemp,Complex(rn/x,0.0));
        cxtemp=Cmul(cxtemp,cxxi);
        cxtemp=Csub(cxtemp,cxxi1);
        cxbn=Cdiv(cxbn,cxtemp);
        /* Augment sums for *qsca and g=<cos(theta)> */
        // *qsca = *qsca + (2.*rn+1.)*(cabs(cxan)**2+cabs(cxbn)**2);
        *qsca = *qsca + (2.*rn+1.)*(Cabs(cxan)*Cabs(cxan)+Cabs(cxbn)*Cabs(cxbn)); 
        *gsca = *gsca + ((2.*rn+1.)/(rn*(rn+1.)))*(cxan.r*cxbn.r+cxan.i*cxbn.i); 
        if (n>1) 
        {
            *gsca = *gsca + ((rn-1.)*(rn+1.)/rn)*(cxan1.r*cxan.r+
            cxan1.i*cxan.i+cxbn1.r*cxbn.r+cxbn1.i*cxbn.i);
        }

        for ( j = 1; j<= nang; j++) 
        {
            jj = 2*nang - j;
            pi[j] = pi1[j];
            tau[j] = rn*amu[j]*pi[j] - (rn+1.E0)*pi0[j];
            p = pow(-1.0,n-1);
            // cxs1[j] = cxs1[j] + fn*(cxan*pi[j]+cxbn*tau[j]);
            cxtemp=Cmul(cxan,Complex(pi[j],0.0));
            cxtemp=Cadd(cxtemp,Cmul(cxbn,Complex(tau[j],0.0)));
            cxtemp=Cmul(Complex(fn,0.0),cxtemp);
            cxs1[j]=Cadd(cxs1[j],cxtemp);
            t = pow(-1.0,n);
            // cxs2[j] = cxs2[j] + fn*(cxan*tau[j]+cxbn*pi[j]);
            cxtemp=Cmul(cxan,Complex(tau[j],0.0));
            cxtemp=Cadd(cxtemp,Cmul(cxbn,Complex(pi[j],0.0)));
            cxtemp=Cmul(Complex(fn,0.0),cxtemp);
            cxs2[j]=Cadd(cxs2[j],cxtemp);

            if (j!=jj) 
            {
                // cxs1[jj] = cxs1[jj] + fn*(cxan*pi(j)*p+cxbn*tau(j)*t);
                cxtemp=Cmul(cxan,Complex(pi[j]*p,0.0));
                cxtemp=Cadd(cxtemp,Cmul(cxbn,Complex(tau[j]*t,0.0)));
                cxtemp=Cmul(Complex(fn,0.0),cxtemp);
                cxs1[jj]=Cadd(cxs1[jj],cxtemp);

                // cxs2[jj] = cxs2[jj] + fn*(cxan*tau(j)*t+cxbn*pi(j)*p);
                cxtemp=Cmul(cxan,Complex(tau[j]*t,0.0));
                cxtemp=Cadd(cxtemp,Cmul(cxbn,Complex(pi[j]*p,0.0)));
                cxtemp=Cmul(Complex(fn,0.0),cxtemp);
                cxs2[jj]=Cadd(cxs2[jj],cxtemp);
            }
        }
        psi0 = psi1;
        psi1 = psi;
        apsi1 = psi1;
        chi0 = chi1;
        chi1 = chi;
        cxxi1 = Complex(apsi1,-chi1);

        /*  For each angle J, compute pi_n+1 from PI = pi_n , PI0 = pi_n-1 */
        for ( j = 1; j<= nang; j++) 
        {
            pi1[j] = ((2.*rn+1.)*amu[j]*pi[j]-(rn+1.)*pi0[j])/rn;
            pi0[j] = pi[j];
        }
    } /*end of big for */

    /*  Have summed sufficient terms. Now compute *qsca,*qext,*qback,and *gsca */
    *gsca = 2.* *gsca/ *qsca;
    *qsca = (2.E0/(x*x))* *qsca;
    *qext = (4.E0/(x*x))*cxs1[1].r;
    *qback = (4.E0/(x*x))*Cabs(cxs1[2*nang-1])*Cabs(cxs1[2*nang-1]);
    free(cxd);
    // *qabs = *qext - *qsca;
    return;
}
void printProgress(double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}
/* 
---------------------------------------------------------
---------------------------------------------------------
If there is an opacity file, read it and get the beta value
---------------------------------------------------------
---------------------------------------------------------
*/
void read_opacity(char *star, Param *param, Sarray *sarray, Warray *warray, RTarray *rtarray)
{
    unsigned int iwav, ig, it, j1, j2;
    double dummy1 = 0., dummy2 = 0., jwav = 0., expterm = 0., jqabs;
    double nu[param->nwav], lstar[param->nwav], lqabs[param->nwav], qbplanck[param->nwav];
    double bplanck_grid[param->nt][param->nwav];
    char filename[1000];

    sprintf(filename, "%s/%s_%.2f_%.2f_%d.dat",star, param->comp, param->smin*1.e4, param->smax*1.e4, param->ng);
    // ------------------------------------
    // I already checked that the file exists. With the proper grain sizes and ng
    // ------------------------------------
    FILE *fichier = fopen(filename, "r" );
    fscanf(fichier,"%d %d", &j1, &j2);
    // ------------------------------------
    // Make some checks on the numbers of grain sizes and wavelengths
    // ------------------------------------
    if (j1 != param->ng)
    {
        printf("The opacity file has the wrong number of grain sizes. Qutting.\n");
        exit(0);
    }
    if (j2 != param->nwav)
    {
        printf("The opacity file has the wrong number of wavelengths. Qutting.\n");
        exit(0);
    }
    for (iwav = 0 ; iwav < param->nwav ; iwav++)
    {
        nu[iwav] = warray[iwav].nu;
        lstar[iwav] = warray[iwav].lstar;
    }
    // ------------------------------------
    // Define the temperature array
    // ------------------------------------
    for (it=0 ; it < param->nt ; it++)
    {
        rtarray[it].btemp_grid = param->tmin * pow((param->tmax/param->tmin),(it*1. / ((param->nt*1.)-1.)));
        for (iwav = 0 ; iwav < param->nwav ; iwav++)
        {
            expterm = exp(warray[iwav].nuhh / rtarray[it].btemp_grid);
            bplanck_grid[it][iwav] = 2.0 * HH * warray[iwav].nu3 / CC2 / (expterm - 1.);
        }
    }
    // ------------------------------------
    // Now, read the file and populate the grain size array as well as compute the beta array
    // ------------------------------------
    for (ig = 0; ig < param->ng ; ig++)
    {
        fscanf(fichier,"%lf", &sarray[ig].gsize);
		sarray[ig].gsize *= 1.e-4;
        for (iwav = 0 ; iwav < param->nwav ; iwav++)
        {
            fscanf(fichier,"%lf %lf", &jwav, &jqabs);
            if (ig == 0) // I make this check only for the first grain size.
            {
                if (fabs(jwav - warray[iwav].wav) > 0.00001)
                {
                    printf("The wavelength grid is different in the opacity file\n");
                    printf("compared to the stellar model. Quitting.\n");
                    exit(0);
                }
            }
			sarray[iwav * param->ng + ig].qabs = jqabs;
            lqabs[iwav] = lstar[iwav] * sarray[iwav * param->ng + ig].qabs;
        }
        // ------------------------------------
        // Now get the radius temperature relationship
        // ------------------------------------
        dummy1 = get_integrate(nu, lqabs, param->nwav) ;
        for (it = 0 ; it < param->nt ; it++)
        {
            for (iwav = 0 ; iwav < param->nwav ; iwav++)
            {
                qbplanck[iwav] = bplanck_grid[it][iwav] * sarray[iwav * param->ng + ig].qabs * M_PI ;
            }
            dummy2 =  get_integrate(nu, qbplanck, param->nwav);
            rtarray[it * param->ng + ig].rt = sqrt(dummy1 / dummy2) * param->dpc * 0.5;
        }
    }
    fclose(fichier);

    /*fichier = fopen("test.dat", "w" );*/
    /*ig = 0;*/
    /*for (it = 0 ; it < param->nt ; it++)*/
    /*{*/
        /*fprintf(fichier, "%lf %lf\n", rtarray[it * param->ng + ig].rt/AU, rtarray[it].btemp_grid);*/
    /*}*/
    /*fclose(fichier);*/

}
/* 
---------------------------------------------------------
---------------------------------------------------------
If there is no opacity file, create it and get the beta array
---------------------------------------------------------
---------------------------------------------------------
*/
void compute_opacity(char *star, Param *param, Sarray *sarray, Warray *warray, RTarray *rtarray)
{
    unsigned int iwav, ig, it, nang_tmp = 2;
    fcomplex cxref, cxs1[mxnang], cxs2[mxnang];
    double x = 0., qext_single = 0., qsca_single = 0., dummy1 = 0., dummy2 = 0., expterm = 0.;
    double nu[param->nwav], lstar[param->nwav], lqabs[param->nwav], qabs[param->nwav], qbplanck[param->nwav];
    float qback = 0., gsca = 0.;
    double bplanck_grid[param->nt][param->nwav];
    char filename[1000];

    sprintf(filename, "%s/%s_%.2f_%.2f_%d.dat",star, param->comp, param->smin * 1.e4, param->smax*1.e4, param->ng);
    FILE *fichier = fopen(filename, "w" );
    fprintf(fichier, "%d %d\n", param->ng, param->nwav);
    // ------------------------------------
    // ------------------------------------
    for (iwav = 0 ; iwav < param->nwav ; iwav++)
    {
        // ------------------------------------
        // Need to make a new array for the integration afterwards
        // ------------------------------------
        nu[iwav] = warray[iwav].nu;
        lstar[iwav] = warray[iwav].lstar;
    }
    // ------------------------------------
    // Define the temperature array
    // ------------------------------------
    for (it=0 ; it < param->nt ; it++)
    {
        rtarray[it].btemp_grid = param->tmin * pow((param->tmax/param->tmin),(it*1. / ((param->nt*1.)-1.)));
        for (iwav = 0 ; iwav < param->nwav ; iwav++)
        {
            expterm = exp(warray[iwav].nuhh / rtarray[it].btemp_grid);
            bplanck_grid[it][iwav] = 2.0 * HH * warray[iwav].nu3 / CC2 / (expterm - 1.);
        }
    }
    // ------------------------------------
    // Define the size array
    // ------------------------------------
    for (ig = 0 ; ig < param->ng; ig++)
    {
        if (param->ng > 1)
        {
            printProgress(1.*ig/param->ng);
            sarray[ig].gsize = param->smin * pow( (param->smax / param->smin), ((1.*ig)/(1. * param->ng - 1.)));
        }
        else
        {
            sarray[ig].gsize = param->smin;            
        }
        fprintf(fichier, "%lf\n", sarray[ig].gsize*1.e4);
        for (iwav = 0 ; iwav < param->nwav ; iwav++)
        {
            // ------------------------------------
            // The optical constant at that given wavelength and the index
            // ------------------------------------
            cxref = Complex(warray[iwav].opt_real, warray[iwav].opt_imag);
            x = 2.e0 * M_PI * sarray[ig].gsize * 1.e4 / warray[iwav].wav ;
            // ------------------------------------
            // Compute Qabs and Qsca. I only need them locally
            // ------------------------------------
            bhmie(x, cxref, nang_tmp, cxs1, cxs2, &qext_single, &qsca_single, &qback, &gsca);
            sarray[iwav * param->ng + ig].qabs = (qext_single - qsca_single);
            lqabs[iwav] = lstar[iwav] * sarray[iwav * param->ng + ig].qabs;
            fprintf(fichier, "%lf %lf\n", warray[iwav].wav, sarray[iwav * param->ng + ig].qabs);
        }
        // ------------------------------------
        // Now get the radius temperature relationship
        // ------------------------------------
        dummy1 = get_integrate(nu, lqabs, param->nwav) ;
        for (it = 0 ; it < param->nt ; it++)
        {
            for (iwav = 0 ; iwav < param->nwav ; iwav++)
            {
                qbplanck[iwav] = bplanck_grid[it][iwav] * qabs[iwav] * M_PI ;
            }
            dummy2 =  get_integrate(nu, qbplanck, param->nwav);
            rtarray[it * param->ng + ig].rt = sqrt(dummy1 / dummy2) * param->dpc * 0.5;
        }
    }
    if (param->ng > 1)
    {
        printf("\n");
    }
    fclose(fichier);        
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Integrate function, trapezoid approximation
---------------------------------------------------------
---------------------------------------------------------
*/
double get_integrate(double *x, double *fx, int len)
{
    double integral = 0.;
    int i = 0;
    for (i = 0 ; i < len-1 ; i++)
    {
        integral = integral + (x[i+1]-x[i]) * (fx[i] + fx[i+1]) * 0.5;
    }
    return integral;
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Get the dust properties
---------------------------------------------------------
---------------------------------------------------------
*/
void get_dust_properties(char *star, Param *param, Sarray *sarray, Warray *warray, RTarray *rtarray, Cla *cla)
{
    unsigned int ig ;
    double dummy = 0., dummy2 = 0., delta_s = 0.;
    // ------------------------------------
    // Read and interpolate properly the optical constant from the comp file
    // ------------------------------------
    get_optical_cst(param, warray);
    // ------------------------------------
    // Prepare the file for the qext file, if it does not exist
    // ------------------------------------
    if (!param->isqfile)
    {
        if (cla->verbose)
        {
            printf("--------------------------------------------------------------------------------\n");
            printf("Computing the opacities ... \n");
        }
        compute_opacity(star, param, sarray, warray, rtarray);
    }
    else
    {
        if (cla->verbose)
        {
            printf("--------------------------------------------------------------------------------\n");
            printf("Reading the opacities ... \n");
        }
        read_opacity(star, param, sarray, warray, rtarray);
    }
    // ------------------------------------
    // Define the ndens array
    // ------------------------------------
    if (param->ng > 1)
    {
        for (ig = 0 ; ig < param->ng; ig++)
        {
            if(ig == 0)
            {
                dummy = sqrt(sarray[ig].gsize * sarray[ig+1].gsize) ;
                delta_s = log10(dummy) - log10(sarray[ig].gsize) ;
            }
            if (ig == param->ng-1)
            {
                dummy = sqrt(sarray[ig-1].gsize * sarray[ig].gsize) ;
                delta_s = log10(sarray[ig].gsize) - log10(dummy) ; 
            }
            if ((ig != 0) && (ig != param->ng-1))
            {
                dummy = sqrt(sarray[ig-1].gsize * sarray[ig].gsize) ; 
                dummy2 = sqrt(sarray[ig].gsize * sarray[ig+1].gsize) ;
                delta_s = log10(dummy2) - log10(dummy) ;
            }
			sarray[ig].ndens = pow(sarray[ig].gsize/sarray[0].gsize, param->grain) * sarray[ig].gsize * delta_s ;
        }
    }
    else
    {
        sarray[0].ndens = 1.e0 ;      
    }
}
/* 
---------------------------------------------------------
---------------------------------------------------------
Interpolation function, in log-space, for x-arrays in decreasing order
---------------------------------------------------------
---------------------------------------------------------
*/
double interpolate_log_down(double *x, double *y, int n, double ix)
{
    int j = binarySearchDown(x, n-1, ix);
    double iy = log(y[j]) + (log(y[j+1]) - log(y[j])) / (log(x[j+1]) - log(x[j])) * (log(ix) - log(x[j]));
    iy = exp(iy);
    return iy;
}
int binarySearchDown(double arr[], int h, double x)
{
    int m, l = 0;
    while (h-l > 1)
    {
        m = l + (h-l)/2;
        if (arr[m] == x)
           return m;
        if (arr[m] > x)
            l = m ;
        else
            h = m;
    }
    return l;
}

