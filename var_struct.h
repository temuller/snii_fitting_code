struct vars_struct {
	// data part
	int *name;
	int *dset;
	int *flt;
	int *ref;
	double *hjd;
	double *val;
	double *err;
	// data pars
	int n_data;
	int n_stars;
	int n_ref;
	int *n_data_per_star;
	int *n_vel_per_star;
	double *min_time;
	double *max_time;
	double *filter_lambda_center;
	double *filter_lambda_width;
	char **sn_name;
	// galaxies pars
	int n_gal;
	int *name_gal;
	bool fit_galaxy;
	char **gal_name;
	// fitting pars
	int n_per_filter;
	int n_per_star;
	int n_phot_bands;
	int n_velo_lines;
	int *fit_star;
	int *fit_vel_star;
	int fit_vel;
	// prior pars
	int n_prior;
	bool* prior_on, *prior_global_on;
	bool prior_on_all;
	int prior_total_on, *prior_subcovar_index;
	double *prior_mean, **prior_covar, **prior_subcovar, **prior_eigenvector, *prior_eigenvalue;
	bool vary_full_reddening_law, vary_only_rv;
	// data blocking
	int n_block_phot_band;
	int *block_phot_band;
	int n_block_star_band;
	int *block_star_band_star,*block_star_band_dset, *block_star_band_band;
	int n_block_star_ref;
	int *block_star_ref_star,*block_star_ref_dset,*block_star_ref_ref;
	int n_block_point;
	int *block_point_star, *block_point_dset, *block_point_band;
	double *block_point_jd;
};
