#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cmpfit-1.2/mpfit.h"
#include "var_struct.h"
#include "read_data.h"
#include "fitujeme.h"
#include "misc.h"
#include "deriv.h"
#include "filter_zeropoints.h"
#include "bolometric.h"

using namespace std;

const double Lsun = 3.9e33;
const double c =2.99792458e10;
const double pc=3.086e18;
const double PI = 3.14159265359;

extern double global_covar[1000][1000];


double mni(vars_struct data, double *a, const int star_no, const int phot_band, const double time_nickel) {
	int star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	double time = a[star_start+0]+time_nickel;
	return bolometric(data,a,star_no,phot_band,time)*Lsun*7.866e-44*exp(( (time-a[star_start+0])  -6.1)/111.26);
}


double log10_mni(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	return log10(mni(data, a, star_no, phot_band, time));
}




void calculate_mni(vars_struct data, double *a, int star_no, bool use_saved_global_covar, char *global_covar_file, double &mni_out, double &l50, double **outcovar, const double time_nickel, const double time_plateau) {
	int n_coef = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*data.n_stars;

	int n_global = data.n_per_filter*data.n_phot_bands+data.n_velo_lines;
	double **covar_global = new double*[n_global];
	for(int i = 0; i < n_global; ++i) {
		covar_global[i] = new double[n_global];
		for (int j=0;j<n_global;j++) covar_global[i][j] = 0.0;
	}

	if (use_saved_global_covar) {read_covar_global(covar_global, n_global, "covar_global.dat");}

	double time, mni_val, mni_sig, *pder_m50=new double[n_coef], *pder_mni = new double[n_coef], m50_val, m50_sig;

	int phot_band = 2;	// dummy phot_band
	int star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;

	time = a[star_start+0]+time_nickel;
	mni_val = log10_mni(data, a, star_no, phot_band, time);
	for (int i=0;i<n_coef;i++) pder_mni[i] = gsl_deriv_central(log10_mni, a, n_coef, i, 1e-4, data, star_no, phot_band, time);

	time = a[star_start+0]+time_plateau;
	m50_val = log10_bolometric(data, a, star_no, phot_band, time);
	for (int i=0;i<n_coef;i++)	pder_m50[i] = gsl_deriv_central(log10_bolometric, a, n_coef, i, 1e-4, data, star_no, phot_band, time);

	for (int i=0; i<2;i++) for (int j=0;j<2;j++) outcovar[i][j] = 0.0;

	if (use_saved_global_covar) {
		for (int i=star_start;i<(star_start+data.n_per_star);i++) for (int j=star_start;j<(star_start+data.n_per_star);j++) {
			outcovar[0][0] += global_covar[i][j]*pder_m50[i]*pder_m50[j];
			outcovar[1][1] += global_covar[i][j]*pder_mni[i]*pder_mni[j];
			outcovar[0][1] += global_covar[i][j]*pder_mni[i]*pder_m50[j];
			outcovar[1][0] += global_covar[i][j]*pder_mni[i]*pder_m50[j];
//			mni_sig += global_covar[i][j]*pder[i]*pder[j];
		}
		for (int i=0;i<n_global;i++) for (int j=0;j<n_global;j++) {
			outcovar[0][0] += covar_global[i][j]*pder_m50[i]*pder_m50[j];
			outcovar[1][1] += covar_global[i][j]*pder_mni[i]*pder_mni[j];
			outcovar[0][1] += covar_global[i][j]*pder_mni[i]*pder_m50[j];
			outcovar[1][0] += covar_global[i][j]*pder_mni[i]*pder_m50[j];
		}
//			mni_sig += covar_global[i][j]*pder[i]*pder[j];
	} else {
		for (int i=0;i<n_coef;i++) for (int j=0;j<n_coef;j++) {
			outcovar[0][0] += global_covar[i][j]*pder_m50[i]*pder_m50[j];
			outcovar[1][1] += global_covar[i][j]*pder_mni[i]*pder_mni[j];
			outcovar[0][1] += global_covar[i][j]*pder_mni[i]*pder_m50[j];
			outcovar[1][0] += global_covar[i][j]*pder_mni[i]*pder_m50[j];
//		mni_sig += global_covar[i][j]*pder[i]*pder[j];
		}
	}


	mni_out = mni_val;
	l50 = m50_val;

	delete[] pder_m50;
	delete[] pder_mni;
	for (int i=0;i<n_global;i++) delete[] covar_global[i];
	delete[] covar_global;
}
