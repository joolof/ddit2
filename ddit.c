#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "fitsio.h"
#include "complex.h"
#include "nrutil.h"
#include "debris.h"
#include "utils.h"
/* 
---------------------------------------------------------
Main program, call all the functions in order
---------------------------------------------------------
*/
// gcc -o ddit debris.c utils.c complex.c nrutil.c -lm -L. -lcfitsio
int main(int argc, char *argv[])
{
	char * star ;
	double seconds ;

	//  ------------------------------------
	// Initialise all the structures
	//  ------------------------------------
	Param param = {0., 0., 0., 0.};
	Cla cla = {0, 0, 0, 0};
	Stellar * stellar ;
}
