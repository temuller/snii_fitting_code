#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "var_struct.h"
#include "math.h"
#include "cmpfit-1.2/mpfit.h"
#include "misc.h"
#include "reddening.h"
#include "fitujeme.h"

extern double global_covar[1000][1000];


void write_covar_snpar(vars_struct data, double *covar, long n_coef) {
	FILE *cov = fopen("covar_snpar.dat", "w");
	int star_start, gal_start;

	int n_index = 11;
	int my_index[11] = { 0, 1, 2, 3, 4, 5,6,7,8,9, 12};

	for (int i=0;i<data.n_stars;i++) {
		star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*i;
		gal_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines;
		if (data.fit_star[i] == 0) {continue;}
		for (int j=0;j<n_index;j++) {
			for (int k=0;k<n_index;k++) fprintf(cov, "%15.10le ", global_covar[star_start+my_index[j]][star_start+my_index[k]]);
			fprintf(cov,"%15.10le ", global_covar[star_start+my_index[j]][gal_start+data.name_gal[i]]);
			fprintf(cov, "\n");
		}
		for (int j=0;j<n_index;j++) fprintf(cov, "%15.10le ", global_covar[gal_start+data.name_gal[i]][star_start+my_index[j]]);
		fprintf(cov, "%15.10le\n", global_covar[gal_start+data.name_gal[i]][gal_start+data.name_gal[i]]);
	}
	fclose(cov);
}


void write_sn_pars_err_full(double *a, vars_struct data, double *uncert, char *filename) {
	FILE *sn_pars = fopen(filename, "w");
	int star_start, nphot, nvel;
	double chi2;
	for (int i=0;i<data.n_stars;i++) {
		if (data.fit_star[i] == 0) {continue;}
		chi2 = calculate_single_chi2(i, a, data, nphot, nvel);
		fprintf(sn_pars, "%-2i %-2i", i, data.name_gal[i]);
		star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+i*data.n_per_star;
		for (int j=0;j<data.n_per_star;j++) fprintf(sn_pars, " %10f %10f", a[star_start+j], uncert[star_start+j]);
		fprintf(sn_pars, " %5i %10f\n", nphot+nvel, chi2);
	}
	fclose(sn_pars);
}



void write_sn_pars_err(double *a, vars_struct data, double *uncert, double *uncert_rescaled) {
	FILE *sn_pars = fopen("sn_pars_err.dat", "w");
	int star_start, gal_start;
	for (int i=0;i<data.n_stars;i++) {
		if (data.fit_star[i] == 0) {continue;}
		fprintf(sn_pars, "%-2i %-2i", i, data.name_gal[i]);
		star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+i*data.n_per_star;
		gal_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines;
		// t0
		fprintf(sn_pars, "   %8.2f %5.2f", a[star_start+0]+data.min_time[i], uncert[star_start+0]);
		// plateau duration
		fprintf(sn_pars, "   %6.2f %5.2f", a[star_start+1], uncert[star_start+1]);
		// plateau width
		fprintf(sn_pars, "   %7.2f %6.2f", a[star_start+2], uncert[star_start+2]);
		// distance modulus
		if (data.fit_galaxy) {
			fprintf(sn_pars, "   %6.3f %5.3f %5.3f", a[gal_start+data.name_gal[i]], uncert[gal_start+data.name_gal[i]], uncert_rescaled[gal_start+data.name_gal[i]]);
		} else {
			fprintf(sn_pars, "   %6.3f %5.3f %5.3f", a[star_start+11], uncert[star_start+11], uncert_rescaled[star_start+11]);
		}
		// reddening
		fprintf(sn_pars, "   %6.3f %5.3f", a[star_start+12], uncert[star_start+12]);
		fprintf(sn_pars, "\n");
	}
	fclose(sn_pars);
}

void write_sn_pars(double *a, vars_struct data) {
	FILE *sn_pars = fopen("sn_pars.dat", "w");
	fprintf(sn_pars, "%i\n", data.n_stars);
	for (int i=0;i<data.n_stars;i++) {
		fprintf(sn_pars, "%i", i);
		fprintf(sn_pars, " %15.10le", a[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+i*data.n_per_star+0]+data.min_time[i]);
		for (int j=1;j<data.n_per_star;j++) fprintf(sn_pars, " %15.10le", a[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+i*data.n_per_star+j]);
		fprintf(sn_pars, "\n");
	}
	fclose(sn_pars);
}



void write_data(double *a, vars_struct data) {

	FILE *phot = fopen("phot.dat", "w");
	FILE *velo = fopen("velo.dat", "w");
	for (int i=0;i<data.n_data; i++) {

		int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*data.name[i];
		if (data.fit_star[data.name[i]] == 0) {continue;}
		if ((data.hjd[i]-a[star_start+0]) > 500) {continue;}

		double theo_mag, theo_vel, wt_exp, rho, temp;
		calculate_mag_vel(data, a, data.name[i], data.hjd[i], data.flt[i], data.flt[i], data.dset[i], theo_mag, theo_vel, wt_exp, rho, temp);
		theo_mag += add_extinction(data, a, data.name[i], data.flt[i])+add_distance(data, a, data.name[i]);

		if (data.dset[i] == 0) {
			// skip if velocity obtained after transition
			//if (data.hjd[i] > (a[star_start+1]+a[star_start+0])) {continue;}
			fprintf(velo, "%f %f %f %i %i %f %f %f\n", data.hjd[i]+data.min_time[data.name[i]], data.val[i], data.err[i], data.flt[i], data.name[i], (data.val[i]-theo_vel)/(wt_exp*data.err[i]), rho, wt_exp);
		}

		if (data.dset[i] == 1) {
			fprintf(phot, "%f %f %f %i %i %f %f %f %f\n", data.hjd[i]+data.min_time[data.name[i]], data.val[i], data.err[i], data.flt[i], data.name[i], (data.val[i]-theo_mag)/data.err[i], rho, temp, wt_exp);
		}
	}
	fclose(velo);
	fclose(phot);
}


void write_model(double *a, vars_struct data) {
	FILE *model_phot = fopen("model_phot.dat", "w");
	FILE *model_velo = fopen("model_velo.dat", "w");
	int n_points = 500, star_start;
	double model_time;
	int flt_list[20] = {0,1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20};
	int n_flt_list = 20;
	fprintf(model_phot, "# Columns: JD");
	for (int i=0;i<n_flt_list;i++) fprintf(model_phot, ", band_%-2i", flt_list[i]);
	fprintf(model_phot, "\n");
	for (int i=0;i<data.n_stars;i++) {
		star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*i;
		for (int j=0;j<n_points;j++) {
			double theo_mag, theo_vel, wt_exp, rho, temp;
			model_time = double(j)/(double(n_points))*550;//(data.max_time[i]-data.min_time[i]) + data.min_time[i];
			fprintf(model_phot, "%15f ", model_time+data.min_time[i]);
			fprintf(model_velo, "%15f ", model_time+data.min_time[i]);

			//for (int k=0;k<data.n_velo_lines;k++) {
				calculate_mag_vel(data, a, i, model_time, 0, 0, 0, theo_mag, theo_vel, wt_exp, rho, temp);
				fprintf(model_velo, "%15f ", theo_vel);
			//}

			for (int k=0;k<n_flt_list;k++) {
				calculate_mag_vel(data, a, i, model_time, flt_list[k], 0, 1, theo_mag, theo_vel, wt_exp, rho, temp);
				theo_mag += add_extinction(data, a, i, flt_list[k])+add_distance(data, a, i);
				fprintf(model_phot, "%15f ", theo_mag);
			}
			fprintf(model_phot, "\n");
			fprintf(model_velo, "\n");
		}

	}
	fclose(model_phot);
	fclose(model_velo);
}
