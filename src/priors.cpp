#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "var_struct.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

void prior_alloc_subcovar(vars_struct data, double **subcovar, double **eigenvector, double *eigenvalue, int *subcovar_index) {
	int n = data.prior_total_on;
	int index[n], count = 0;
	gsl_matrix * gsl_subcovar = gsl_matrix_alloc(n,n);
  	gsl_vector *gsl_eval = gsl_vector_alloc (n);
	gsl_matrix *gsl_evec = gsl_matrix_alloc (n, n);
	gsl_eigen_symmv_workspace *gsl_w = gsl_eigen_symmv_alloc(n);

	printf("n: %i\n", n);

	// determine which indices are switched on
	for (int i=0;i<data.n_per_star;i++) {
		if (data.prior_on[i]) {
			index[count] = i;
			subcovar_index[count] = index[count];
			count++;
		}
	}
	for (int i=0;i<n;i++) {
		for (int j=0;j<n;j++) {
			subcovar[i][j] = data.prior_covar[index[i]][index[j]];
			gsl_matrix_set(gsl_subcovar, i, j, subcovar[i][j]);
//			printf("%15.10le ", subcovar[i][j]);
		}
//		printf("\n");
	}

	gsl_eigen_symmv(gsl_subcovar, gsl_eval, gsl_evec, gsl_w);
	for (int i=0;i<n;i++) {
		eigenvalue[i] = gsl_vector_get(gsl_eval, i);
		for (int j=0;j<n;j++) eigenvector[i][j] = gsl_matrix_get(gsl_evec, i, j);
	}
	gsl_vector_free(gsl_eval);
	gsl_matrix_free(gsl_evec);
	gsl_matrix_free(gsl_subcovar);
	gsl_eigen_symmv_free(gsl_w);
}


void calculate_prior(double *a, vars_struct &data) {
	int count[data.n_per_star];
	for (int i=0;i<data.n_per_star;i++) {
		data.prior_mean[i] = 0.0;
		for (int j=0;j<data.n_per_star;j++) data.prior_covar[i][j] = 0.0;
		count[i] = 0;
	}
	int star_start;
	for (int i=0; i<data.n_stars;i++) {
		if (data.fit_star[i] == 0) {continue;}
		if (data.n_vel_per_star[i] < 5) {continue;}
		star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*i;
/*
		// no prior on explosion time
		// plateau and transition width are just from all stars
		for (int j=1;j<=2;j++) {
			par_mean[j] += a[star_start +j];
			count[j]++;
		}
		// mean of velocity params only from stars that have velocities
		if (data.n_vel_per_star[i] > 0)
		for (int j=3; j<=5;j++) {
			par_mean[j] += a[star_start +j];
			count[j]++;
		}
		// the rest from all stars

		for (int j=6;j<=10;j++) {
			par_mean[j] += a[star_start +j];
			count[j]++;
		}
	*/
		for (int j=0;j<data.n_per_star;j++) {
			data.prior_mean[j] += a[star_start+j];
			count[j]++;
		}
	}

	// calulcate the mean
	for (int i=0; i<data.n_per_star;i++) {
		if (count[i] == 0) {continue;}
		data.prior_mean[i] = data.prior_mean[i]/double(count[i]);
	}

	// now calculate the covariance matrix
	for (int i=0;i<data.n_per_star;i++) count[i] = 0;
	for (int i=0; i<data.n_stars;i++) {
		if (data.fit_star[i] == 0) {continue;}
		if (data.n_vel_per_star[i] < 5) {continue;}
		count[0]++;
		star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*i;
		for (int j=0;j<data.n_per_star;j++)	{
			for (int k=0; k<data.n_per_star;k++) data.prior_covar[j][k] += (a[star_start+j]-data.prior_mean[j])*(a[star_start+k]-data.prior_mean[k]);
		}
	}
	for (int i=0;i<data.n_per_star;i++) for (int j=0;j<data.n_per_star;j++) {
		data.prior_covar[i][j] = data.prior_covar[i][j]/double(count[0]-1);
	}


}


void write_prior_star(vars_struct data) {
	FILE *prior = fopen("sn_prior.dat", "w");
	fprintf(prior, "%i\n", data.n_per_star);
	for (int i=0;i<data.n_per_star;i++) fprintf(prior, "%i %15.10le\n", i, data.prior_mean[i]);
	for (int i=0;i<data.n_per_star;i++) {
		for (int j=0;j<data.n_per_star;j++) fprintf(prior, "%15.10le ", data.prior_covar[i][j]);
		fprintf(prior, "\n");
	}
	fclose(prior);
}


void read_prior_star(double *prior_mean, double **prior_covar, vars_struct data, const char *filein) {
	FILE *dov = fopen(filein, "r");
	char dummy[1000];
	int n_per_star;
	fgets(dummy, 1000, dov);
	sscanf(dummy, "%i", &n_per_star);
	if (n_per_star != data.n_per_star) {
		fprintf(stderr, "read_prior_star: number of pars per star does not match the data\n");
		exit(0);
	}
	double a0;
	int t;
	for (int i=0; i<n_per_star;i++) {
		fgets(dummy, 1000, dov);
		sscanf(dummy, "%i %le", &t, &a0);
		prior_mean[t] = a0;
	}
	for (int i=0; i<n_per_star;i++) {
		for (int j=0;j<n_per_star;j++) {
			fscanf(dov, "%le", &a0);
			prior_covar[i][j] = a0;
		}
	}
}
