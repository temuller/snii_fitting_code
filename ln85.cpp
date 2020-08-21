#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "var_struct.h"
#include "fitujeme.h"
#include "deriv.h"
#include "bolometric.h"
#include "read_data.h"
#include "mni.h"
#include "physical_parameters.h"


const double Lsun = 3.9e33;
const double c =2.99792458e10;
const double pc=3.086e18;
const double PI = 3.14159265359;

extern double global_covar[1000][1000];



void calculate_ln85(vars_struct data, double *a, const int star_no) {
		int n_coef = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*data.n_stars;
		int n_global = data.n_per_filter*data.n_phot_bands+data.n_velo_lines;

	double **covar = new double*[n_global];
	for (int i=0;i<n_global;i++) covar[i] = new double[n_global];
	read_covar_global(covar, n_global, "covar_global.dat");


		int star_start;

		double time_nickel =200;
		double time_plateau = 50;

		double en, *vals = new double[5], *times = new double[5];

		double **pder = new double*[5], **outcovar = new double*[5];
		for (int i=0;i<5;i++) {
			pder[i] = new double[n_coef];
			outcovar[i] = new double[5];
		}


		int phot_band = 2;
		FILE *ven = fopen("ln85_err_rescale.dat", "w");


			star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
			int dist_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines + data.name_gal[star_no];
			vals[0] = log10_mni(data, a, star_no, phot_band, time_nickel);
			vals[1] = log10_bolometric_plateau(data, a, star_no, phot_band, time_plateau);
			vals[2] = ln85_e(data, a, star_no, phot_band, 0);
			vals[3] = ln85_m(data, a, star_no, phot_band, 0);
			vals[4] = ln85_r(data, a, star_no, phot_band, 0);

			for (int i=0;i<n_coef;i++) {
				pder[0][i] = gsl_deriv_central(log10_mni, a, n_coef, i, 1e-4, data, star_no, phot_band, time_nickel);
				pder[1][i] = gsl_deriv_central(log10_bolometric_plateau, a, n_coef, i, 1e-4, data, star_no, phot_band, time_plateau);
				pder[2][i] = gsl_deriv_central(ln85_e, a, n_coef, i, 1e-4, data, star_no, phot_band, 0);
				pder[3][i] = gsl_deriv_central(ln85_m, a, n_coef, i, 1e-4, data, star_no, phot_band, 0);
				pder[4][i] = gsl_deriv_central(ln85_r, a, n_coef, i, 1e-4, data, star_no, phot_band, 0);
			}

			for (int k=0;k<5;k++) for (int l=0;l<5;l++) outcovar[k][l] = 0.0;
			for (int k=0;k<5;k++) for (int l=0;l<5;l++) {
				// take care of uncertainties due to fit, comes from the fit covar. matrix global_covar
				for (int i=star_start;i<(star_start+data.n_per_star);i++) for (int j=star_start;j<(star_start+data.n_per_star);j++) outcovar[k][l] += global_covar[i][j]*pder[k][i]*pder[l][j];
				// now add covariances from the global parameters (should be tiny) using saved matrix in covar
				for (int i=0;i<n_global;i++) for (int j=0;j<n_global;j++) outcovar[k][l] += covar[i][j]*pder[k][i]*pder[l][j];
//				for (int j=0;j<n_coef;j++) outcovar[k][l] += covar[i][j]*pder[k][i]*pder[l][j];
			}
			fprintf(ven, "%i %f %f ", star_no, a[star_start+0], data.max_time[star_no]-data.min_time[star_no]-a[star_start+0]);
			for (int i=0;i<5;i++) fprintf(ven, "%f ", vals[i]);
			for (int i=0;i<5;i++) for (int j=0;j<5;j++) fprintf(ven, "%le ", outcovar[i][j]);
			fprintf(ven, "\n");

		fclose(ven);


		for (int i=0;i<n_global;i++) delete[] covar[i];
		delete[] covar;
		delete[] pder;
}

void calculate_popov93(vars_struct data, double *a, const int star_no) {
		int n_coef = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*data.n_stars;
		int n_global = data.n_per_filter*data.n_phot_bands+data.n_velo_lines;

	double **covar = new double*[n_global];
	for (int i=0;i<n_global;i++) covar[i] = new double[n_global];
	read_covar_global(covar, n_global, "covar_global.dat");


		int star_start;

		double time_nickel =200;
		double time_plateau = 50;

		double en, *vals = new double[5], *times = new double[5];

		double **pder = new double*[5], **outcovar = new double*[5];
		for (int i=0;i<5;i++) {
			pder[i] = new double[n_coef];
			outcovar[i] = new double[5];
		}


		int phot_band = 2;
		FILE *ven = fopen("popov93_err_rescale.dat", "w");


			star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
			int dist_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines + data.name_gal[star_no];
			vals[0] = log10_mni(data, a, star_no, phot_band, time_nickel);
			vals[1] = log10_bolometric_plateau(data, a, star_no, phot_band, time_plateau);
			vals[2] = popov93_e(data, a, star_no, phot_band, 0);
			vals[3] = popov93_m(data, a, star_no, phot_band, 0);
			vals[4] = popov93_r(data, a, star_no, phot_band, 0);

			for (int i=0;i<n_coef;i++) {
				pder[0][i] = gsl_deriv_central(log10_mni, a, n_coef, i, 1e-4, data, star_no, phot_band, time_nickel);
				pder[1][i] = gsl_deriv_central(log10_bolometric_plateau, a, n_coef, i, 1e-4, data, star_no, phot_band, time_plateau);
				pder[2][i] = gsl_deriv_central(popov93_e, a, n_coef, i, 1e-4, data, star_no, phot_band, 0);
				pder[3][i] = gsl_deriv_central(popov93_m, a, n_coef, i, 1e-4, data, star_no, phot_band, 0);
				pder[4][i] = gsl_deriv_central(popov93_r, a, n_coef, i, 1e-4, data, star_no, phot_band, 0);
			}

			for (int k=0;k<5;k++) for (int l=0;l<5;l++) outcovar[k][l] = 0.0;
			for (int k=0;k<5;k++) for (int l=0;l<5;l++) {
				// take care of uncertainties due to fit, comes from the fit covar. matrix global_covar
				for (int i=star_start;i<(star_start+data.n_per_star);i++) for (int j=star_start;j<(star_start+data.n_per_star);j++) outcovar[k][l] += global_covar[i][j]*pder[k][i]*pder[l][j];
				// now add covariances from the global parameters (should be tiny) using saved matrix in covar
				for (int i=0;i<n_global;i++) for (int j=0;j<n_global;j++) outcovar[k][l] += covar[i][j]*pder[k][i]*pder[l][j];
//				for (int j=0;j<n_coef;j++) outcovar[k][l] += covar[i][j]*pder[k][i]*pder[l][j];
			}
			fprintf(ven, "%i %f %f ", star_no, a[star_start+0], data.max_time[star_no]-data.min_time[star_no]-a[star_start+0]);
			for (int i=0;i<5;i++) fprintf(ven, "%f ", vals[i]);
			for (int i=0;i<5;i++) for (int j=0;j<5;j++) fprintf(ven, "%le ", outcovar[i][j]);
			fprintf(ven, "\n");

		fclose(ven);


		for (int i=0;i<n_global;i++) delete[] covar[i];
		delete[] covar;
		delete[] pder;
}

