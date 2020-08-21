#include <math.h>
#include "var_struct.h"
#include "fitujeme.h"
#include "bolometric.h"
#include "mni.h"

double ln85_e(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	double theo_mag, theo_vel, wt_exp, rho, temp, t_midp;
	int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	t_midp = a[star_start+0]+0.5*a[star_start+1];
	calculate_mag_vel(data, a, star_no, t_midp, 2, 0, 1, theo_mag, theo_vel, wt_exp, rho, temp);
	theo_mag -= 30;
	//printf("%i %f %f %f\n", star_no, theo_mag, a[star_start+0], a[star_start+1]);
	return 0.135*theo_mag + 2.34*log10(a[star_start+1]) + 3.13*log10(theo_vel/1000.0) - 3.205;
}

double ln85_m(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	double theo_mag, theo_vel, wt_exp, rho, temp, t_midp;
	int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	t_midp = a[star_start+0]+0.5*a[star_start+1];
	calculate_mag_vel(data, a, star_no, t_midp, 2, 0, 1, theo_mag, theo_vel, wt_exp, rho, temp);
	theo_mag -= 30;
	//printf("%i %f %f %f\n", star_no, theo_mag, a[star_start+0], a[star_start+1]);
	return 0.234*theo_mag + 2.91*log10(a[star_start+1]) + 1.96*log10(theo_vel/1000.0) - 1.829;
}

double ln85_r(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	double theo_mag, theo_vel, wt_exp, rho, temp, t_midp;
	int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	t_midp = a[star_start+0]+0.5*a[star_start+1];
	calculate_mag_vel(data, a, star_no, t_midp, 2, 0, 1, theo_mag, theo_vel, wt_exp, rho, temp);
	theo_mag -= 30;
	//printf("%i %f %f %f\n", star_no, theo_mag, a[star_start+0], a[star_start+1]);
	return -0.572*theo_mag - 1.07*log10(a[star_start+1]) -2.74*log10(theo_vel/1000.0) - 3.350;
}


double popov93_e(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	double theo_mag, theo_vel, wt_exp, rho, temp, t_midp;
	int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	t_midp = a[star_start+0]+0.5*a[star_start+1];
	calculate_mag_vel(data, a, star_no, t_midp, 2, 0, 1, theo_mag, theo_vel, wt_exp, rho, temp);
	theo_mag -= 30;
	//printf("%i %f %f %f\n", star_no, theo_mag, a[star_start+0], a[star_start+1]);
	return 0.4*theo_mag + 4.0*log10(a[star_start+1]) + 5.0*log10(theo_vel/1000.0) - 4.311 + 1.0;
}

double popov93_m(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	double theo_mag, theo_vel, wt_exp, rho, temp, t_midp;
	int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	t_midp = a[star_start+0]+0.5*a[star_start+1];
	calculate_mag_vel(data, a, star_no, t_midp, 2, 0, 1, theo_mag, theo_vel, wt_exp, rho, temp);
	theo_mag -= 30;
	//printf("%i %f %f %f\n", star_no, theo_mag, a[star_start+0], a[star_start+1]);
	return 0.4*theo_mag + 4.0*log10(a[star_start+1]) + 3.0*log10(theo_vel/1000.0) - 2.089;
}

double popov93_r(vars_struct data, double *a, const int star_no, const int phot_band, const double time) {
	double theo_mag, theo_vel, wt_exp, rho, temp, t_midp;
	int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*star_no;
	t_midp = a[star_start+0]+0.5*a[star_start+1];
	calculate_mag_vel(data, a, star_no, t_midp, 2, 0, 1, theo_mag, theo_vel, wt_exp, rho, temp);
	theo_mag -= 30;
	//printf("%i %f %f %f\n", star_no, theo_mag, a[star_start+0], a[star_start+1]);
	return -0.8*theo_mag - 2.0*log10(a[star_start+1]) -4.0*log10(theo_vel/1000.0) - 4.278;
}

