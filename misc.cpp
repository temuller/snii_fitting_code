#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cmpfit-1.2/mpfit.h"
#include "var_struct.h"
#include "fitujeme.h"
#include "reddening.h"
#include "misc.h"


double calculate_single_chi2(int fit_star_no, double *a, vars_struct data, int &single_nphot, int &single_nvel) {
	double dy[data.n_data];
	double chi2 = 0;
	bool blocked;
	single_nphot = 0;
	single_nvel = 0;
		int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*fit_star_no;
		double wt_exp;
		for (int i=0;i<data.n_data;i++) {
			if (data.name[i] != fit_star_no) {continue;}
			dy[i] = 0;

			if (data.fit_star[data.name[i]] == 0) {continue;}
			if ((data.hjd[i]-a[star_start+0]) < 0) {continue;}
			if ((data.hjd[i]-a[star_start+0]) > 500) {continue;}

		double theo_mag, theo_vel, rho, temp;

		calculate_mag_vel(data, a, data.name[i], data.hjd[i], data.flt[i], data.flt[i], data.dset[i], theo_mag, theo_vel, wt_exp, rho, temp);

		if (data.dset[i] == 0) {
			// only Fe II so far
			if (data.flt[i] != 0) {continue;}
			// do not fit if velocity fitting shut off
			if (data.fit_vel == 0) {continue;}
			// do not fit if velocity fitting shut off for this particular star
			if (data.fit_vel_star[data.name[i]] == 0) {continue;}

			dy[i] = (data.val[i] - theo_vel)/(data.err[i]/wt_exp);
			single_nvel++;

		}

		if (data.dset[i] == 1) {
			theo_mag += add_extinction(data, a, data.name[i], data.flt[i]);
			theo_mag += add_distance(data, a, data.name[i]);

			dy[i] = (data.val[i]-theo_mag)/data.err[i];
			single_nphot++;
		}
		chi2 += dy[i]*dy[i];
	}
	return chi2;
}




void set_cardelli(double *a, vars_struct data, double rv) {
	for (int i=0;i<data.n_phot_bands;i++) {
		//printf("%i %f %f\n", i, data.filter_lambda_center[i], reddening_coef(data.filter_lambda_center[i], rv));
		a[data.n_per_filter*i + 0] = reddening_coef(data.filter_lambda_center[i], rv);
	}
}


void shutoff_extinction_coef(mp_par *pars, vars_struct data) {
	for (int i=0;i<data.n_phot_bands;i++) pars[data.n_per_filter*i+0].fixed = 1;
}

void vary_full_reddening_law(mp_par *pars, vars_struct &data) {
	data.vary_full_reddening_law = true;
	data.vary_only_rv = false;
	for (int i=0;i<data.n_phot_bands;i++) pars[data.n_per_filter*i+0].fixed = 0;
}

void vary_only_rv(mp_par *pars, vars_struct &data) {
	data.vary_full_reddening_law = false;
	data.vary_only_rv = true;
	for (int i=0;i<data.n_phot_bands;i++) pars[data.n_per_filter*i+0].fixed = 1;
	pars[data.n_per_filter*2+0].fixed = 0;
	pars[data.n_per_filter*2+0].limited[0] = 1; pars[data.n_per_filter*2+0].limits[0] = 1.0;
	pars[data.n_per_filter*2+0].limited[1] = 1; pars[data.n_per_filter*2+0].limits[1] = 7.0;
}


void set_star_coef (double *a, vars_struct data, int coef_no, double val) {
	// sets value of a parameter to the same for all stars
	for (int i=0;i<data.n_stars;i++) {
		a[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+i*data.n_per_star+coef_no] = val;
	}
}

void set_filter_coef (double *a, vars_struct data, int coef_no, double val) {
	// sets value of a parameter to the same for all stars
	for (int i=0;i<data.n_phot_bands;i++) {
		a[data.n_per_filter*i+coef_no] = val;
	}
}


void shutoff_star_coef (mp_par *pars, vars_struct data, int coef_no) {
	// do not fit one parameter for all stars
	for (int i=0;i<data.n_stars;i++) {
		pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+i*data.n_per_star+coef_no].fixed = 1;
	}
}



void shutoff_vel_star(int star_no, mp_par *pars, vars_struct data) {
	data.fit_vel_star[star_no] = 0;
}

void shutoff_vel_all(mp_par *pars, vars_struct data) {
	for (int i=0;i<data.n_velo_lines;i++) pars[data.n_per_filter*data.n_phot_bands+i].fixed = 1;
	data.fit_vel = 0;
}

void fit_only_star(int star_no, mp_par *pars, vars_struct data) {
	data.fit_star[star_no] = 1;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+0].fixed = 0;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+1].fixed = 0;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+2].fixed = 0;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+3].fixed = 0;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+4].fixed = 0;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+5].fixed = 0;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+6].fixed = 0;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+7].fixed = 0;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+8].fixed = 0;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+9].fixed = 0;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+10].fixed = 1;	// do not fit temperature during exponential decay
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+11].fixed = 0;
	pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+12].fixed = 0;
}


void shutoff_star(int star_no, mp_par *pars, vars_struct data) {
	for (int i=0; i<data.n_per_star;i++) {
		pars[data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+star_no*data.n_per_star+i].fixed = 1;
	}
	data.fit_star[star_no] = 0;
}


void shutoff_vel_except_iron  (mp_par *pars, vars_struct data) {
	// fit only Fe II velocities (index 0)
	for (int i=1; i<data.n_velo_lines; i++) {
		pars[data.n_per_filter*data.n_phot_bands+i].fixed = 1;
	}
}

void shutoff_sn_pars (mp_par *pars, vars_struct data) {
	// do not fit any supernova-specific parameters

	for (int i=0; i<data.n_stars; i++) {
		int cur_index = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+i*data.n_per_star;
		for (int j=0;j<data.n_per_star;j++) pars[cur_index+j].fixed = 1;
	}

}


void shutoff_global (mp_par *pars, vars_struct data) {
	// do not fit global parameters

	for (int i=0; i<data.n_phot_bands; i++) for (int j=0;j<data.n_per_filter;j++) {
		pars[data.n_per_filter*i+j].fixed = 1;
	}

	for (int i=0; i<data.n_velo_lines; i++) {
		pars[data.n_per_filter*data.n_phot_bands+i].fixed = 1;
	}
}

void shutoff_temp (mp_par *pars, vars_struct data) {
	// do not fit temperature coefficients

	for (int i=0; i<data.n_phot_bands; i++) {
		pars[data.n_per_filter*i+2].fixed = 1;
		pars[data.n_per_filter*i+3].fixed = 1;
		pars[data.n_per_filter*i+4].fixed = 1;
	}
}


void shutoff_temp2 (mp_par *pars, vars_struct data) {
	// do not fit higher temperature coefficients than linear

	for (int i=0; i<data.n_phot_bands; i++) {
		pars[data.n_per_filter*i+3].fixed = 1;
		pars[data.n_per_filter*i+4].fixed = 1;
	}
}

void shutoff_temp3 (mp_par *pars, vars_struct data) {
	// do not fit higher temperature coefficients than quadratic
	for (int i=0; i<data.n_phot_bands; i++) {
		pars[data.n_per_filter*i+4].fixed = 1;
	}
}

void initialize_pars(mp_par *pars, double *a, vars_struct data, int n_coef) {
	printf("init pars\n");
	for (int i=0;i<n_coef;i++) {
		pars[i].fixed = 1;
		a[i] = 0;
	}

	int n_data, n_stars, n_phot_bands, n_velo_lines, n_gal, n_per_filter, n_per_star;
	int *dset, *flt;
	n_data = data.n_data;
	n_stars = data.n_stars;
	n_phot_bands = data.n_phot_bands;
	n_velo_lines = data.n_velo_lines;
	n_gal  = data.n_gal;
	n_per_filter = data.n_per_filter;
	n_per_star = data.n_per_star;
	dset = data.dset;
	flt = data.flt;

	printf("n_data: %i, n_star: %i\n", n_data, n_stars);

	int run_index = 0;
	// first do the temperature coefficients
	for (int i=0; i<n_phot_bands; i++) {
		a[n_per_filter*i + 0] = 3.1;
		// ratio of total to selective extinction R_i = A_i/E(B-V)
		pars[n_per_filter*i + 0].limited[0] = 1; pars[n_per_filter*i + 0].limits[0] = 0.0;
		pars[n_per_filter*i + 0].limited[1] = 1; pars[n_per_filter*i + 0].limits[1] = 15.0;
		// absolute coeff;
		a[n_per_filter*i + 1] = 31.0;
		// linear coefficients;
		a[n_per_filter*i + 2] = 5.0;
		// second order (whatever it is)
		a[n_per_filter*i + 3] = 0.0;
		// thrid order (whatever it is)
		a[n_per_filter*i + 4] = 0.0;
	}
	run_index += n_per_filter*n_phot_bands;

	// add velocity projection factors;
	for (int i=0; i<n_velo_lines;i++) {
		a[run_index + i] = 1.0;
		pars[run_index + i].limited[0] = 1;
		pars[run_index + i].limits[0] = 0.0;
	}
	run_index += n_velo_lines;

	// galaxy distances
	for (int i=0;i<n_gal;i++) {
		a[run_index +i] = 0.0;
		pars[run_index+i].limited[0] = 1; pars[run_index+i].limits[0] = -10;
		pars[run_index+i].limited[1] = 1; pars[run_index+i].limits[1] = 10;
	}
	run_index += n_gal;


	for (int i=0;i<n_stars;i++) {
		// current index
		int cur_ind = run_index + i*n_per_star;
		// basic time information
		a[cur_ind + 0] = data.min_time[i]-2;	// JD0
		pars[cur_ind + 0].fixed = 0;
		pars[cur_ind + 0].limited[1] = 1; pars[cur_ind + 0].limits[1] = -0.01; //data.min_time[i];
		pars[cur_ind + 0].limited[0] = 1; pars[cur_ind + 0].limits[0] = -200; //data.min_time[i];
		a[cur_ind + 1] = 130;	// Transition time in days
		pars[cur_ind + 1].fixed = 0;
		pars[cur_ind + 1].limited[0] = 1; pars[cur_ind + 1].limits[0] = 0.0;
		pars[cur_ind + 1].limited[1] = 1; pars[cur_ind + 1].limits[1] = 140.0;
		a[cur_ind + 2] = 5;	// Transition width in days
		pars[cur_ind + 2].fixed = 0;
		pars[cur_ind + 2].limited[0] = 1; pars[cur_ind + 2].limits[0] = 0.1;
		pars[cur_ind + 2].limited[1] = 1; pars[cur_ind + 2].limits[1] = 50.0;

		//vel_theo =  a[0]*(JD-JD0)^a[1]+a[2]
		a[cur_ind + 3] = 3.7e4;
		pars[cur_ind + 3].fixed = 0;
		pars[cur_ind + 3].limited[0] = 1; pars[cur_ind +3].limits[0] = 10.0;
		a[cur_ind + 4] = -0.58;
		pars[cur_ind + 4].fixed = 0;
		pars[cur_ind + 4].limited[0] = 1; pars[cur_ind +4].limits[0] = -1.0;
		pars[cur_ind + 4].limited[1] = 1; pars[cur_ind +4].limits[1] = 0.0;
		a[cur_ind + 5] = 504;
		pars[cur_ind + 5].fixed = 0;
		pars[cur_ind +5].limited[0] = 1; pars[cur_ind + 5].limits[0] = 0.0;

		// linear radius decay
		a[cur_ind + 6] = -4e-3;
		pars[cur_ind + 6].fixed = 0;
		pars[cur_ind + 6].limited[1] = 1; pars[cur_ind + 6].limits[1] = -0.001;
		a[cur_ind + 7] = 10;
		pars[cur_ind + 7].fixed = 0;

		// temperature during expansion
		a[cur_ind + 8] = -0.005;
		pars[cur_ind + 8].fixed = 0;
		pars[cur_ind + 8].limited[1] = 1; pars[cur_ind + 8].limits[1] = 0.0;
		a[cur_ind + 9] = 0.1;
		pars[cur_ind + 9].fixed = 0;

		// temperature during decay
		a[cur_ind + 10] = -0.4;
		pars[cur_ind + 10].fixed = 1;

		// relative distance modulus
		a[cur_ind + 11] = 0.0;
		pars[cur_ind + 11].fixed = 0;

		// E(B-V)
		a[cur_ind + 12] = 0.03;
		pars[cur_ind + 12].fixed = 0;
		pars[cur_ind + 12].limited[0] = 1; pars[cur_ind + 12].limits[0] = 0.0;
	}

	// shut off everything where there is no enough data
	for (int i=0;i<n_data; i++) {
		if (dset[i] == 0) {
			pars[n_per_filter*n_phot_bands + flt[i]].fixed = 0;
		}
		if ((dset[i] == 1) && (data.fit_star[data.name[i]] == 1)) {
			// switch coeffs;
			pars[n_per_filter*flt[i]+0].fixed = 0;
			pars[n_per_filter*flt[i]+1].fixed = 0;
			pars[n_per_filter*flt[i]+2].fixed = 0;
			pars[n_per_filter*flt[i]+3].fixed = 0;
			pars[n_per_filter*flt[i]+4].fixed = 0;
		}

	}
	// velocity relative to Fe II
	pars[n_per_filter*n_phot_bands+0].fixed = 1;
}


double perform_fit(vars_struct data, mp_par *pars, double *a, double *uncert, double *covar, int n_coef, int code) {
// does the fit and returns results. Set "code" for different verbosity of the ouput


		mp_result result;

		memset(&result, 0, sizeof(result));
		mp_config config;
		memset(&config, 0, sizeof(config));
		config.maxiter=500;
		result.xerror = uncert;
		result.covar = covar;

//for (long i=0;i<n_coef*n_coef;i++) printf("%i %15.10le\n", i, covar[i]);

		int status;
//		printf("zde %i %i %f\n", data.n_data, n_coef, a[6]);
		status = mpfit(ridicka, data.n_data+data.n_prior, n_coef, a, pars, &config, (void *) &data, &result);
//		printf("zde\n");
		if (code >= 1) 	printf("*** testlinfit status = %d, niter = %i\n", status, result.niter);
		if (code >= 2)	{
			printf("orig CHI2    = %f    (%d DOF)\n", result.orignorm, result.nfunc-result.nfree);
			printf("  CHI-SQUARE = %f    (%d DOF)\n", result.bestnorm, result.nfunc-result.nfree);
		}
		if (code >= 3) {
			printf("        NPAR = %d\n", result.npar);
			printf("       NFREE = %d\n", result.nfree);
			printf("     NPEGGED = %d\n", result.npegged);
			printf("     NITER = %d\n", result.niter);
			printf("      NFEV = %d\n", result.nfev);
		}
		if (code >= 4) {
			printf("\n");
	   	for (int i=0; i<result.npar; i++) {
		      printf("  P[%d] = %f +/- %f\n", i, a[i], result.xerror[i]);
    	}
	}
	return result.bestnorm;

}


