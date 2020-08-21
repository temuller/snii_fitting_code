// reads in data. This is the only part of the code in C++ and not C.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "var_struct.h"
#include "read_data.h"
#include "priors.h"

using namespace std;


void read_covar_global(double **covar, int n_coef, char *filename) {
	FILE *dov = fopen(filename, "r");
	double dummy;
	for (int i=0;i<n_coef;i++) {
		for (int j=0;j<n_coef;j++) {
			fscanf(dov, "%le", &dummy);
			covar[i][j] = dummy;
		}
	}
	fclose(dov);
}


void read_covar_full(double **covar, int n_coef, char *filename) {
	FILE *dov = fopen(filename, "r");
	double dummy;
	for (int i=0;i<n_coef;i++) {
		for (int j=0;j<n_coef;j++) {
			fscanf(dov, "%le", &dummy);
			covar[i][j] = dummy;
		}
	}
	fclose(dov);
}


void initialize_data (vars_struct &data, char *filename) {
		const int maxdata = 15000;
		int *name = new int[maxdata], *dset = new int[maxdata], *flt = new int[maxdata], *ref = new int[maxdata];
		for (int i=0;i<maxdata;i++) name[i] = 0;
		double *hjd = new double[maxdata], *val = new double[maxdata], *err = new double[maxdata];
		int n_data, n_stars=1, n_phot_bands=28, n_velo_lines=15, n_ref;
		const int n_per_star = 13, n_per_filter=5;
		read_data(ref, dset, flt, hjd, val, err, n_data, n_phot_bands, n_velo_lines, n_ref, filename);
		printf("ndata: %i nref %i\n", n_data, n_ref);

		double *filter_lambda_center = new double[n_phot_bands], *filter_lambda_width = new double[n_phot_bands], *prior_mean = new double[n_per_star];
		double **prior_covar = new double*[n_per_star];
		for(int i = 0; i < n_per_star; ++i) prior_covar[i] = new double[n_per_star];

		bool *prior_on = new bool[n_per_star], *prior_global_on = new bool[n_per_filter];

		for (int i=0;i<n_phot_bands;i++) {
			filter_lambda_center[i] = 0;
			filter_lambda_width[i] = 0;
		}

		int *fit_star = new int[n_stars];
		for (int i=0;i<n_stars; i++) fit_star[i] = 1;

		double *min_time = new double[n_stars], *max_time = new double[n_stars];
		for (int i=0; i<n_stars; i++) {
			min_time[i] = 99999;
			max_time[i] = 0;
		}
		for (int i=0; i<n_data; i++) {
			if (min_time[name[i]] > hjd[i]) {min_time[name[i]] = hjd[i];}
			if (max_time[name[i]] < hjd[i]) {max_time[name[i]] = hjd[i];}
		}
		// reduce the HJD to the first observed point for numerical reasons
		for (int i=0;i<n_data; i++) {
			hjd[i] = hjd[i] - min_time[name[i]];
		}

		int *fit_vel_star = new int[n_stars];
		for (int i=0;i<n_stars;i++) {fit_vel_star[i] = 1;}
		int fit_vel = 1;
/*
		int *n_data_per_star = new int[n_stars], *n_vel_per_star = new int[n_stars];
		for (int i=0;i<n_stars;i++) {n_data_per_star[i] = 0; n_vel_per_star[i] = 0;}
		for (int i=0;i<n_data;i++) {
			n_data_per_star[name[i]]++;
			if ((dset[i] == 0) && (flt[i] == 0)) n_vel_per_star[name[i]]++;
		}

		int *n_data_per_filter = new int[n_phot_bands];
		for (int i=0;i<n_phot_bands;i++) n_data_per_filter[i] = 0;
		for (int i=0;i<n_data;i++) {
			if (dset[i] == 1) n_data_per_filter[flt[i]]++;
		}
*/
//		for (int i=0;i<n_stars;i++) {printf("%i %i\n", i, n_data_per_star[i]);}


		printf("%i %i %i %i %i\n", n_data, n_stars, n_ref, n_phot_bands, n_velo_lines);
//		for (int i =0; i<n_data; i++) {
//			if (err[i] <= 0) printf("XXXXXXXXXXXXXXXXXXXXXX\n");
//			printf("%i %3i %3i %3i %10.4f %8.4f %7.4f\n", i, name[i], dset[i], flt[i], hjd[i], val[i], err[i]);
//		}



	data.name = name;
	data.ref = ref;
	data.dset = dset;
	data.flt = flt;
	data.hjd = hjd;
	data.val = val;
	data.err = err;
	data.n_data = n_data;
	data.n_stars = n_stars;
	data.n_phot_bands = n_phot_bands;
	data.n_velo_lines = n_velo_lines;
	data.n_ref = n_ref;
	data.fit_star = fit_star;
	data.min_time = min_time;
	data.max_time = max_time;
	data.fit_vel_star = fit_vel_star;
	data.fit_vel = fit_vel;
	data.n_per_star = n_per_star;
	data.n_per_filter = n_per_filter;
	data.filter_lambda_center = filter_lambda_center;
	data.filter_lambda_width = filter_lambda_width;
//	data.n_data_per_star = n_data_per_star;
//	data.n_vel_per_star = n_vel_per_star;
	data.vary_full_reddening_law = false;
	data.vary_only_rv = false;
	data.prior_on = prior_on;
	data.prior_global_on = prior_global_on;

	read_prior_star(prior_mean, prior_covar, data, "sn_prior.dat");
	data.prior_mean = prior_mean;
	data.prior_covar = prior_covar;

}



int process_line (const string &line, int &ref, int &data, int &flt, double &hjd, double &val, double &err) {
		stringstream ss;
		int stat;
		ss << line;
		ss >> ref >> data >> flt >> hjd >> val >> err;
		stat = 1-int(ss.fail());
		return stat;
}

int read_block_point(int *star, int *dset, int *band, double *jd, const char *filein) {
	FILE *dov = fopen(filein, "r");
	char dummy[1000];
	int n_blocked, a0, a1,a2, n=0;
	double a3;
	fgets(dummy, 1000, dov);
	sscanf(dummy, "%i", &n_blocked);
	for (int i=0; i<n_blocked; i++) {
		fgets(dummy, 1000, dov);
		sscanf(dummy, "%i %i %i %lf", &a0, &a1, &a2, &a3);
		star[n] = a0;
		dset[n] = a1;
		band[n] = a2;
		jd[n] = a3;
		n++;
	}
	return n_blocked;
}


int read_block_star_ref(int *star, int *dset, int *ref, const char *filein) {
	FILE *dov = fopen(filein, "r");
	char dummy[1000];
	int n_blocked, a0, a1,a2, n=0;
	fgets(dummy, 1000, dov);
	sscanf(dummy, "%i", &n_blocked);
	for (int i=0; i<n_blocked; i++) {
		fgets(dummy, 1000, dov);
		sscanf(dummy, "%i %i %i", &a0, &a1, &a2);
		star[n] = a0;
		dset[n] = a1;
		ref[n] = a2;
		n++;
	}
	return n_blocked;
}


int read_block_star_band(int *star, int *dset, int *band, const char *filein) {
	FILE *dov = fopen(filein, "r");
	char dummy[1000];
	int n_blocked, a0, a1,a2, n=0;
	fgets(dummy, 1000, dov);
	sscanf(dummy, "%i", &n_blocked);
	for (int i=0; i<n_blocked; i++) {
		fgets(dummy, 1000, dov);
		sscanf(dummy, "%i %i %i", &a0, &a1, &a2);
		star[n] = a0;
		dset[n] = a1;
		band[n] = a2;
		n++;
	}
	return n_blocked;
}


void read_sn_pars(double *a, vars_struct data, const char *filein) {
	FILE *dov = fopen(filein, "r");
	char dummy[1000];
	int n_stars;
	fgets(dummy, 1000, dov);
	sscanf(dummy, "%i", &n_stars);
	if (n_stars != data.n_stars) {
		fprintf(stderr, "read_sn_pars: number of stars does not match the data\n");
		exit(0);
	}
	double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12;
	int t, cur_ind;
	for (int i=0; i<n_stars;i++) {
		fgets(dummy, 1000, dov);
		sscanf(dummy, "%i %le %le %le %le %le %le %le %le %le %le %le %le %le", &t, &a0, &a1, &a2, &a3, &a4, &a5, &a6, &a7, &a8, &a9, &a10, &a11, &a12);
		cur_ind = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+i*data.n_per_star;
		a[cur_ind+0] = a0;
		a[cur_ind+1] = a1;
		a[cur_ind+2] = a2;
		a[cur_ind+3] = a3;
		a[cur_ind+4] = a4;
		a[cur_ind+5] = a5;
		a[cur_ind+6] = a6;
		a[cur_ind+7] = a7;
		a[cur_ind+8] = a8;
		a[cur_ind+9] = a9;
		a[cur_ind+10] = a10;
		a[cur_ind+11] = a11;
		a[cur_ind+12] = a12;
	}
}

void read_global(double *a,  vars_struct data, const char *filein) {
	FILE *dov = fopen(filein, "r");
	char dummy[200];
	int n_phot_bands, n_velo_lines, n_gal;
	fgets(dummy, 200, dov);
	sscanf(dummy, "%i %i %i", &n_phot_bands, &n_velo_lines, &n_gal);
	double d1, d2, d3, d4, d5, lam, dlam;
	int di;
	for (int i=0;i<n_phot_bands; i++) {
		fgets(dummy, 200, dov);
		sscanf(dummy, "%i %le %le %le %le %le %le %le", &di, &lam, &dlam, &d1, &d2, &d3, &d4, &d5);
		data.filter_lambda_center[i] = lam;
		data.filter_lambda_width[i] = dlam;
		a[data.n_per_filter*i+0] = d1;
		a[data.n_per_filter*i+1] = d2;
		a[data.n_per_filter*i+2] = d3;
		a[data.n_per_filter*i+3] = d4;
		a[data.n_per_filter*i+4] = d5;
	}
	for (int i=0;i<n_velo_lines; i++) {
		fgets(dummy, 200, dov);
		sscanf(dummy, "%i %le", &di, &d1);
		a[data.n_per_filter*n_phot_bands+i] = d1;
	}
	fclose(dov);
}


void read_data (int *ref, int *dset, int *flt, double *hjd, double *val, double *err, int &n_data, int &n_phot_bands, int &n_velo_lines, int &n_ref, const char *filein) {
	int t_n_data =0, line_stat;
	int t_data, t_flt, t_ref;
	double t_hjd, t_val, t_err;
	string line;
	ifstream file_input;
	file_input.open(filein);
	while ( file_input.good() ){
		getline (file_input,line);
		line_stat = process_line(line, t_ref, t_data, t_flt, t_hjd, t_val, t_err);
		if (line_stat == 1) {
			//printf("%i %i %f %f %f\n", t_data, t_flt, t_hjd, t_val,t_err);
			ref[t_n_data] = t_ref;
			dset[t_n_data] = t_data;
			flt[t_n_data] = t_flt;
			hjd[t_n_data] = t_hjd;
			val[t_n_data] = t_val;
			err[t_n_data] = t_err;
			t_n_data++;
			if (n_ref < (t_ref+1)) {n_ref = t_ref+1;}
		}
	}
	n_data = t_n_data;
	file_input.close();
	file_input.clear();

}
