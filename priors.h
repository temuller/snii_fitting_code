void calculate_prior(double *a, vars_struct &data);
void write_prior_star(vars_struct data);
void read_prior_star(double *prior_mean, double **prior_covar, vars_struct data, const char *filein);
void prior_alloc_subcovar(vars_struct data, double **subcovar, double **eigenvector, double *eigenvalue, int *subcovar_index);