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
#include "reddening.h"
#include "priors.h"
#include "deriv.h"
#include "filter_zeropoints.h"
#include "bolometric.h"
#include "mni.h"

using namespace std;

const double Lsun = 3.9e33;
const double c =2.99792458e10;
const double pc=3.086e18;
const double PI = 3.14159265359;

extern double global_covar[1000][1000];

double change_t_nick(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	double time_plateau = 50;
	int star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	return log10_mni(data, a, star_no, phot_band, a[star_start+0]+time) - log10_bolometric(data, a, star_no, phot_band, a[star_start+0]+time_plateau);
}

double compute_tau(vars_struct data, double *a, const int star_no, const double time) {
	int star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	double wt_exp = 1.0/(1.0+exp((time-a[star_start+0]-a[star_start+1])/a[star_start+2]));
	double wt_dec = 1.0-wt_exp;
	return (a[star_start+8]*(time-a[star_start+0])+a[star_start+9])*wt_exp+  a[star_start+10]*wt_dec;
}

double t_tau0(vars_struct data, double *a, const int star_no) {
	double guess_t[2], old_tau[2], guess_tau, guess_time, tol=1e-4;
	int star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	guess_t[0] = a[star_start+0]-100;
	guess_t[1] = a[star_start+0]+200;
	old_tau[0] = compute_tau(data,a,star_no,guess_t[0]);
	old_tau[1] = compute_tau(data,a,star_no,guess_t[1]);

	while ( fabs(guess_t[1]-guess_t[0]) > tol) {
		guess_time = 0.5*(guess_t[0]+guess_t[1]);
		guess_tau = compute_tau(data,a,star_no,guess_time);
		if (guess_tau < 0) {guess_t[1] = guess_time;} else {guess_t[0] = guess_time;}
	}
	return guess_t[0];
}


double mni_tau0(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	double time_nickel = 200;
	int star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	return log10_mni(data, a, star_no, phot_band, a[star_start+0]+time_nickel) - log10_bolometric(data, a, star_no, phot_band, a[star_start+0]+t_tau0(data,a,star_no));
}


double change_t_plat(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	double time_nickel = 200;
	int star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	return log10_mni(data, a, star_no, phot_band, a[star_start+0]+time_nickel) - log10_bolometric(data, a, star_no, phot_band, a[star_start+0]+time);
}



void mni_dependence(vars_struct data, double *a, const int star_no) {
	int n_coef = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*data.n_stars;
	int n_global = data.n_per_filter*data.n_phot_bands+data.n_velo_lines;

	double **covar = new double*[n_global];
	for (int i=0;i<n_global;i++) covar[i] = new double[n_global];
	read_covar_global(covar, n_global, "covar_global.dat");


	int star_start;

	int n_nick = 6, n_plat=9;

	double time_nickel[n_nick];
	time_nickel[0] = 150;
	time_nickel[1] = 200;
	time_nickel[2] = 250;
	time_nickel[3] = 300;
	time_nickel[4] = 350;
	time_nickel[5] = 400;

	double time_plateau[n_plat];
	time_plateau[0] = 10;
	time_plateau[1] = 20;
	time_plateau[2] = 30;
	time_plateau[3] = 40;
	time_plateau[4] = 50;
	time_plateau[5] = 60;
	time_plateau[6] = 70;
	time_plateau[7] = 80;
	time_plateau[8] = 90;


	double *vals = new double[n_nick];
	double **pder = new double*[n_nick], **outcovar = new double*[n_nick];
	for (int i=0;i<n_nick;i++) {
		pder[i] = new double[n_coef];
		outcovar[i] = new double[n_nick];
	}


	int phot_band = 2;
	FILE *ven = fopen("mni_nickel.dat", "w");

		star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
		for (int i=0;i<n_nick;i++) {
			vals[i] = change_t_nick(data, a, star_no, phot_band, time_nickel[i]);
			for (int j=0;j<n_coef;j++)	pder[i][j] = gsl_deriv_central(change_t_nick, a, n_coef, j, 1e-4, data, star_no, phot_band, time_nickel[i]);
		}

		for (int k=0;k<n_nick;k++) for (int l=0;l<n_nick;l++) outcovar[k][l] = 0.0;
		for (int k=0;k<n_nick;k++) for (int l=0;l<n_nick;l++) {
				// take care of uncertainties due to fit, comes from the fit covar. matrix global_covar
				for (int i=star_start;i<(star_start+data.n_per_star);i++) for (int j=star_start;j<(star_start+data.n_per_star);j++) outcovar[k][l] += global_covar[i][j]*pder[k][i]*pder[l][j];
				// now add covariances from the global parameters (should be tiny) using saved matrix in covar
				for (int i=0;i<n_global;i++) for (int j=0;j<n_global;j++) outcovar[k][l] += covar[i][j]*pder[k][i]*pder[l][j];
			//		for (int i=0;i<n_coef;i++) for (int j=0;j<n_coef;j++) outcovar[k][l] += covar[i][j]*pder[k][i]*pder[l][j];
		}

		fprintf(ven, "%i %f %f ", star_no, a[star_start+0], data.max_time[star_no]-data.min_time[star_no]-a[star_start+0]);
		for (int i=0;i<n_nick;i++) fprintf(ven, "%f %f ", vals[i], sqrt(outcovar[i][i]));
		fprintf(ven, "\n");
	fclose(ven);

	delete[] vals;
	for (int i=0;i<n_nick;i++) {
		delete[] pder[i];
		delete[] outcovar[i];
	}
	delete[] pder;
	delete[] outcovar;




	vals = new double[n_plat];
	pder = new double*[n_plat];
	outcovar = new double*[n_plat];
	for (int i=0;i<n_plat;i++) {
		pder[i] = new double[n_coef];
		outcovar[i] = new double[n_plat];
	}

	ven = fopen("mni_plateau.dat", "w");

		star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
		for (int i=0;i<n_plat;i++) {
			vals[i] = change_t_plat(data, a, star_no, phot_band, time_plateau[i]);
			for (int j=0;j<n_coef;j++)	pder[i][j] = gsl_deriv_central(change_t_plat, a, n_coef, j, 1e-4, data, star_no, phot_band, time_plateau[i]);
		}

		for (int k=0;k<n_plat;k++) for (int l=0;l<n_plat;l++) outcovar[k][l] = 0.0;
		for (int k=0;k<n_plat;k++) for (int l=0;l<n_plat;l++) {
						// take care of uncertainties due to fit, comes from the fit covar. matrix global_covar
				for (int i=star_start;i<(star_start+data.n_per_star);i++) for (int j=star_start;j<(star_start+data.n_per_star);j++) outcovar[k][l] += global_covar[i][j]*pder[k][i]*pder[l][j];
						// now add covariances from the global parameters (should be tiny) using saved matrix in covar
				for (int i=0;i<n_global;i++) for (int j=0;j<n_global;j++) outcovar[k][l] += covar[i][j]*pder[k][i]*pder[l][j];
		}
	//	for (int j=0;j<n_coef;j++) outcovar[k][l] += covar[i][j]*pder[k][i]*pder[l][j];

		fprintf(ven, "%i %f %f ", star_no, a[star_start+0], data.max_time[star_no]-data.min_time[star_no]-a[star_start+0]);
		for (int i=0;i<n_plat;i++) fprintf(ven, "%f %f ", vals[i], sqrt(outcovar[i][i]));
		fprintf(ven, "\n");
	fclose(ven);

	ven = fopen("mni_tau0.dat", "w");
		star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;

		vals[0] = mni_tau0(data, a, star_no, phot_band, 0);
		for (int j=0;j<n_coef;j++)	pder[0][j] = gsl_deriv_central(mni_tau0, a, n_coef, j, 1e-4, data, star_no, phot_band, 0);


		outcovar[0][0] = 0.0;
//		for (int i=0;i<n_coef;i++) for (int j=0;j<n_coef;j++)
				// take care of uncertainties due to fit, comes from the fit covar. matrix global_covar
				for (int i=star_start;i<(star_start+data.n_per_star);i++) for (int j=star_start;j<(star_start+data.n_per_star);j++) outcovar[0][0] += global_covar[i][j]*pder[0][i]*pder[0][j];
				// now add covariances from the global parameters (should be tiny) using saved matrix in covar
				for (int i=0;i<n_global;i++) for (int j=0;j<n_global;j++) outcovar[0][0] += covar[i][j]*pder[0][i]*pder[0][j];

//		outcovar[0][0] += covar[i][j]*pder[0][i]*pder[0][j];

		fprintf(ven, "%i %f %f %f %f %f", star_no, a[star_start+0], data.max_time[star_no]-data.min_time[star_no]-a[star_start+0], t_tau0(data,a,star_no), vals[0], sqrt(outcovar[0][0]));
		fprintf(ven, "\n");
	fclose(ven);

}
