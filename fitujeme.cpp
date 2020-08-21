#include <math.h>
#include <stdio.h>
#include "var_struct.h"
#include "reddening.h"
#include "cmpfit-1.2/mpfit.h"
#include "misc.h"

double add_extinction(vars_struct data, double *a, const int star_no, const int phot_band) {
	int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	int filt_start = data.n_per_filter*phot_band;
	double use_ri;
	// normally use Cardelli law
	use_ri = reddening_coef(data.filter_lambda_center[phot_band], 3.3);
	// varying only RV within Cardelli law
	if (data.vary_only_rv) {
		use_ri = reddening_coef(data.filter_lambda_center[phot_band], a[data.n_per_filter*2+0]);
	}
	// varying full reddening law, maintaing R_B-R_V=1
	if (data.vary_full_reddening_law) {
		use_ri = a[filt_start+0];
		if (phot_band == 1) use_ri = a[data.n_per_filter*2+0]+1.0;
	}

	return a[star_start+12]*use_ri;
}

double add_distance(vars_struct data, double *a, const int star_no) {
	int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	int gal_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines;
	double mu;
	// add distance modulus
	if (data.fit_galaxy) {
		mu = a[gal_start+data.name_gal[star_no]];
	} else {
		mu = a[star_start+11];
	}
	return mu;
}


void calculate_mag_vel(vars_struct data, double *a, const int star_no, const double time, const int phot_band, const int vel_band, const int dset, double &mag, double &velo, double &wt_exp, double &rho, double &temp) {
		int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
		double wt_dec, velo_theo, r_exp, model;


		// alternative from Olivares, Hamuy et al.
		wt_exp = 1.0/(1.0+exp((time-a[star_start+0]-a[star_start+1])/a[star_start+2]));
		wt_dec = 1.0-wt_exp;
		velo_theo = a[star_start+3]*pow(time-a[star_start+0], a[star_start+4])+a[star_start+5];
		velo = a[data.n_per_filter*data.n_phot_bands+vel_band]*velo_theo;

		if (dset == 0) {
			mag = 0.0;
			temp = 0.0;
			return;
		}

		int filt_start = data.n_per_filter*phot_band;
		r_exp = velo_theo*(time-a[star_start+0]);
		rho = log10(wt_exp*r_exp*r_exp + wt_dec*pow(10.0, a[star_start+6]*(time-a[star_start+0]-a[star_start+1]) + a[star_start+7]));

		temp = (a[star_start+8]*(time-a[star_start+0])+a[star_start+9])*wt_exp+  a[star_start+10]*wt_dec;
		model = 30.0+ a[filt_start+1] - 2.5*rho - 2.5*a[filt_start+2]*temp - 2.5*a[filt_start+3]*temp*temp/2.0- 2.5*a[filt_start+4]*temp*temp*temp/6.0;

		mag = model;
		return;
}



// this is the main fitting routine
long ridicka(int n_tot, int n_coef, double *a, double *dy, double **derivs, void *vars) {

	struct vars_struct *data = (struct vars_struct *) vars;
	double *hjd, *val, *err,  *prior_mean;
	int *name, *dset, *flt, *fit_star, n_phot_bands, n_velo_lines, n_stars, n_gal, fit_vel, *fit_vel_star, n_per_star, n_per_filter, n_data, n_prior, *prior_subcovar_index;
	bool prior_on_all;
	// read data from the structure
	name = data->name;
	dset = data->dset;
	flt = data->flt;
	hjd = data->hjd;
	val = data->val;
	err = data->err;
	n_stars = data->n_stars;
	n_phot_bands = data->n_phot_bands;
	n_velo_lines = data->n_velo_lines;
	n_gal = data->n_gal;
	fit_star = data->fit_star;
	fit_vel = data->fit_vel;
	fit_vel_star = data->fit_vel_star;
	n_per_star = data->n_per_star;
	n_per_filter = data->n_per_filter;
	n_data = data->n_data;
	n_prior = data->n_prior;
	prior_mean = data->prior_mean;
	prior_subcovar_index = data->prior_subcovar_index;
	prior_on_all = data->prior_on_all;


	for (int i=0;i<n_data; i++) {

		//printf("pracuju %i\n", i);

		int star_start = n_per_filter*n_phot_bands+n_velo_lines+n_gal+n_per_star*0;
		double wt_exp;
		dy[i] = 0;

		if ((hjd[i]-a[star_start+0]) < 0) {continue;}
		if ((hjd[i]-a[star_start+0]) > 500) {continue;}

		double theo_mag, theo_vel, rho, temp;

		calculate_mag_vel(*data, a, 0, hjd[i], flt[i], flt[i], dset[i], theo_mag, theo_vel, wt_exp, rho, temp);

		if (dset[i] == 0) {
			// only Fe II so far
			if (flt[i] != 0) {continue;}
			// do not fit if velocity fitting shut off
			if (fit_vel == 0) {continue;}
			// do not fit if velocity fitting shut off for this particular star
			if (fit_vel_star[name[i]] == 0) {continue;}

			dy[i] = (val[i] - theo_vel)/(err[i]/wt_exp);

		}

		if (dset[i] == 1) {
			theo_mag += add_extinction(*data, a, 0, flt[i]);
			theo_mag += add_distance(*data, a, 0);

			dy[i] = (val[i]-theo_mag)/err[i];
		}
//		printf("%i %f\n", i, dy[i]);

	}


	// priors

	for (int i=0;i<n_prior;i++) dy[n_data+i] = 0;

	if (prior_on_all) {
		int filt_start;
		for (int i=0;i<n_phot_bands;i++) {
			filt_start = n_per_filter*i;
			if (data->vary_full_reddening_law) dy[n_data + filt_start + 0] = (a[filt_start + 0] - reddening_coef(data->filter_lambda_center[i], 3.3))/(0.1*reddening_coef(data->filter_lambda_center[i], 3.3));
		}

	int star_start;
		star_start = n_per_filter*n_phot_bands+n_velo_lines+n_gal+n_per_star*0;
		for (int j=0;j< data->prior_total_on;j++) for (int k=0;k< data->prior_total_on;k++) {
				dy[n_data+star_start+j] += data->prior_eigenvector[k][j]*(a[star_start+prior_subcovar_index[k]]-prior_mean[prior_subcovar_index[k]])/sqrt(data->prior_eigenvalue[j]);
		}



	}

	return 0;


}
