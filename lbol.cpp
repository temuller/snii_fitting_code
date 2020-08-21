#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cmpfit-1.2/mpfit.h"
#include "var_struct.h"
#include "fitujeme.h"
#include "read_data.h"
#include "deriv.h"
#include "filter_zeropoints.h"
#include "bolometric.h"

using namespace std;

const double Lsun = 3.9e33;
const double c =2.99792458e10;
const double pc=3.086e18;
const double PI = 3.14159265359;

extern double global_covar[1000][1000];



void calculate_lbol (vars_struct data, double *a) {
	int n_coef = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*data.n_stars;

	int n_global = data.n_per_filter*data.n_phot_bands+data.n_velo_lines;
	double **covar_global = new double*[n_global];
	for(int i = 0; i < n_global; ++i) {
		covar_global[i] = new double[n_global];
		for (int j=0;j<n_global;j++) covar_global[i][j] = 0.0;
	}

	read_covar_global(covar_global, n_global, "covar_global.dat");

	FILE *ven = fopen("lbol_err.dat", "w");
	double time, mag, sigma_snpars, sigma_global, *pder=new double[n_coef];

	int star_no = 0;

	int phot_band = 2;
	int star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;

	for (int kk = 0; kk<400;kk++) {
		time = kk;
		mag = bolometric(data, a, star_no, phot_band, time);
		for (int i=0;i<n_coef;i++) {
			pder[i] = gsl_deriv_central(bolometric, a, n_coef, i, 1e-4, data, star_no, phot_band, time);
		}

		sigma_snpars = 0.0;
		sigma_global = 0.0;
		for (int i=star_start;i<n_coef;i++) for (int j=star_start;j<n_coef;j++) sigma_snpars += global_covar[i][j]*pder[i]*pder[j];
		for (int i=0;i<n_global;i++) for (int j=0;j<n_global;j++) sigma_global += covar_global[i][j]*pder[i]*pder[j];

		fprintf(ven, "%f %le %le %le %le\n", time+data.min_time[0], mag, sqrt(sigma_snpars), sqrt(sigma_global), sqrt(sigma_snpars+sigma_global));
	}
	fclose(ven);
	delete[] covar_global;
}
