int ridicka(int n_data, int n_coef, double *a, double *dy, double **derivs, void *vars);
double add_extinction(vars_struct data, double *a, const int star_no, const int phot_band);
double add_distance(vars_struct data, double *a, const int star_no);
void calculate_mag_vel(vars_struct data, double *a, const int star_no, const double time, const int phot_band, const int vel_band, const int dset, double &mag, double &velo, double &wt_exp, double &rho, double &temp);
