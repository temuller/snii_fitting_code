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

const double Lsun = 3.9e33;
const double c =2.99792458e10;
const double pc=3.086e18;
const double PI = 3.14159265359;

double return_magnitude(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	double theo_mag, theo_vel, wt_exp, rho, temp;
	calculate_mag_vel(data, a, star_no, time,phot_band, 0, 1, theo_mag, theo_vel, wt_exp, rho, temp);
	return theo_mag;
}

double bolometric(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	// phot_band is dummy, used for compatibility
	// filters entering bolometric calculation in the order of increasing central wavelength
	const int nbol_flt = 12;
	int bol_flt[nbol_flt] = {13, 14, 15, 16, 0, 1, 2, 3, 4, 5, 6, 7};
	double *flt_zero = new double[data.n_phot_bands], *fluxes = new double[nbol_flt];
	// define filter zeropoints
	define_filter_zeropoints(flt_zero, data.filter_lambda_center, data.n_phot_bands);
	double theo_mag, theo_vel, wt_exp, rho, temp;
	for (int i=0;i<nbol_flt;i++) {
		calculate_mag_vel(data, a, star_no, time, bol_flt[i], 0, 1, theo_mag, theo_vel, wt_exp, rho, temp);
		fluxes[i] = flt_zero[bol_flt[i]]*pow(10.0, -0.4*(theo_mag-30.0))*1e4; //factor 1e4 is due to conversion from Angstrom to microns
	}
	int i_bolo_start=0;
	//if (temp < -0.06) {i_bolo_start = 4;}	// include Swift only if tau > -0.05 (where it is well-defined)
	double bolo = 0.0;
	for (int i=i_bolo_start;i<nbol_flt-1;i++) {
		bolo += 0.5*(fluxes[i+1]+fluxes[i])*(data.filter_lambda_center[bol_flt[i+1]]-data.filter_lambda_center[bol_flt[i]]);
	}
	bolo += fluxes[nbol_flt-1]*data.filter_lambda_center[bol_flt[nbol_flt-1]]/3.0;
	bolo *= 4.0*PI*(10.0*pc)*(10.0*pc)/Lsun;
	delete[] flt_zero;
	delete[] fluxes;
	return (bolo);
}

double log10_bolometric(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	int star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	return log10(bolometric(data, a, star_no, phot_band, time));	
}


double log10_bolometric_plateau(vars_struct data, double *a, const int star_no, const int phot_band, const double time_pl) {
	int star_start =data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	double time = a[star_start+0]+time_pl;
	return log10_bolometric(data, a, star_no, phot_band, time);
	
}



