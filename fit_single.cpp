#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cmpfit-1.2/mpfit.h"
#include "var_struct.h"
#include "read_data.h"
#include "write_data.h"
#include "fitujeme.h"
#include "misc.h"
#include "reddening.h"
#include "lbol.h"
#include "mni.h"
#include "priors.h"
#include "ln85.h"
#include "mni_dependence.h"
#include "physical_parameters.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


using namespace std;

#include <omp.h>

double global_covar[1000][1000];

int main (int argc, char *argv[]) {
	struct vars_struct data;
	initialize_data(data, argv[1]);
	for (int i=0;i<data.n_data;i++) {
		if (data.err[i]< 0.02) data.err[i] = 0.02;
	}

	data.n_phot_bands = 28;
	data.n_velo_lines = 15;
	data.n_gal = 0;

	int n_coef = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*data.n_stars;
	data.n_prior = n_coef;
	printf("ncoef: %i\n", n_coef);

	mp_par pars[n_coef];
	double a[n_coef], orig_a[n_coef], a_rescale[n_coef];
	memset(&pars[0], 0, sizeof(pars));
	initialize_pars(pars, a, data, n_coef);

	read_global(a, data, "global.dat");
	for (int i=0;i<n_coef;i++) orig_a[i] = a[i];

//	for (int i=0;i<n_coef;i++) printf("%i %15.5le %i %i %i %15.5le %15.5le\n", i, a[i], pars[i].fixed, pars[i].limited[0], pars[i].limited[1], pars[i].limits[0], pars[i].limits[1]);

	int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal;


	// SN name is used in the output files	
	string sn_name;
	cout << "Name of the SN: ";
	cin  >> sn_name;

        // create directory with the SN name
    	mkdir(sn_name.c_str(), 0777); 

	

// Here initial values on parameters and their constraints are determined

// Parameter 0: t_exp - setting limits here is of importance to get reasonable results with little data
//	Uncomment this line to specify the upper limit on t_exp (ie the latest time the star could have exploded
//	The time in the code is measured in days relative to the JD of the first observation of each star (to prevent numerical problems with big values of JD
//	The sample lines below sets the upper limit on t_exp to 56533.875
//
//pars[star_start+0].limited[1] = 1;
//pars[star_start+0].limits[1]=56530.0944922-data.min_time[0];

//pars[star_start+0].limited[0] = 1;
//pars[star_start+0].limits[0]=56528.0944922-data.min_time[0];

//
//  Similarly, a lower limit is set by a uncommenting:
//
//		pars[star_start+0].limited[1] = 1;
//		pars[star_start+0].limits[0]=56503.0-data.min_time[0];
//
//	If you know exactly the explosion time, you can set it fixed to a specific value by uncommenting these two lines
//
//  	pars[star_start+0].fixed = 1;
//		a[star_start+0] = 342-data.min_time[0];
//
//	Otherwise, t_exp is set to 8 days before the upper limit:
        
        float texp_val, texp_err;
	printf("Please enter t_exp (0 if no constrain wanted): ");
    	scanf("%f", &texp_val);
    	if(texp_val!=0){
		printf("Please enter std: ");
    		scanf("%f", &texp_err);
		pars[star_start+0].limited[1] = 1;
		pars[star_start+0].limits[1]= (texp_val-data.min_time[0]) + texp_err;
		pars[star_start+0].limited[0] = 1;
		pars[star_start+0].limits[0]= (texp_val-data.min_time[0]) - texp_err;
		a[star_start+0] = (texp_val-data.min_time[0]);
	}
        else{a[star_start+0] = -10;}
        //else{ a[star_start+0] = pars[star_start+0].limits[1]-1; }
        //a[star_start+0] = -10;

//
//	Note that in all cases the initial value of t_exp in a[] needs to be in between the pars.limits, ie. this must be satisfied:
//		pars[star_start+0].limits[0] < a[star_start+0] < pars[star_start].limits[1]
//	Similar inequalities apply to other parameters.
//
//	Some of these parameters already have some general limits on their values set in function initialize_pars in file misc.cpp
//
// Parameter 1: t_P duration of the plateau. You can set upper and lower limits just like for t_exp, just make sure the starting value is between these two limits otherwise fitting won't work
        a[star_start+1] = 80;
	/*float tp_val, tp_err;
	printf("Please enter t_plateau (0 if no limit wanted): ");
    	scanf("%f", &tp_val);
    	if(tp_val!=0){
		printf("Please enter std: ");
    		scanf("%f", &tp_err);
		pars[star_start+1].limited[1] = 1;
		pars[star_start+1].limits[1]= tp_val+tp_err;
		pars[star_start+1].limited[0] = 1;
		pars[star_start+1].limits[0]= tp_val-tp_err;
		a[star_start+1] = tp_val;
	}
        else{a[star_start+1] = 80;}*/

// Parameter 2: t_w duration of the transition width
	a[star_start + 2] = 5.7203710127e+00;
// Parameter 3: omega_0 absolute scale of the velocity
	a[star_start + 3] = 80000;
// Parameter 4: omega_1 velocity exponent
	a[star_start + 4] = -0.7;
// Parameter 5: omega_2 overall velocity offset, usually remains 0
	a[star_start + 5] =  0.0;
// Parameter 6: gamma_0 slope of the exponential decay
	a[star_start + 6] =  -4.3116491519e-03;
// Parameter 7: gamma_1 normalization of the exponential decay
	a[star_start + 7] = 1.0946588903e+01;
//	pars[star_start + 7].limited[1] = 1;
//	pars[star_start + 7].limits[1] = 12;
// Parameter 8: alpha_0 slope of the temperature decrease
	a[star_start + 8] = -5.4171738803e-03;
//	pars[star_start+8].fixed=1;
// Parameter 9: alpha_1 temperature at t_exp
	a[star_start + 9] = 1.1807361799e-01;
// Parameter 10: alpha_thin temperature during exponential decay, should remain fixed at -0.4 (by default)
	a[star_start + 10]= -4.0000000000e-01;
// Parameter 11: mu distance modulus with 30.00 substracted, ie. distance modulus of 34.5 is entered as 4.5. Can also be fixed, free by default
	a[star_start + 11] = 0.0;
//	pars[star_start + 11].fixed=1;
//	pars[star_start+11].fixed =1;
// Parameter 12: E(B-V)
	a[star_start + 12] = 0.6;



	shutoff_vel_except_iron(pars, data);

	data.prior_on_all = true;
	for (int i=0;i<data.n_per_star;i++) data.prior_on[i] = false;

	// switch on prior on plateau duration
	data.prior_on[1] = true;
	// switch on prior on transition duration
	data.prior_on[2] = true;
	// switch on velocity prior
	for (int i=3;i<=5;i++) data.prior_on[i] = true;
	// switch on prior for linear radius decay
	data.prior_on[6] = true;
	data.prior_on[7] = true;
	// switch on prior for temperature changes
	data.prior_on[8] = true;
	data.prior_on[9] = true;


	int prior_total_on = 0;
	for (int i=0;i<data.n_per_star;i++) if (data.prior_on[i]) {prior_total_on++;}
	data.prior_total_on = prior_total_on;

	double **subcovar = new double*[prior_total_on];
	double **eigenvector = new double*[prior_total_on];
	int prior_subcovar_index[prior_total_on];
	for(int i = 0; i<prior_total_on; ++i) {
		//printf("%i %i\n", i, prior_total_on);
		subcovar[i] = new double[prior_total_on];
		eigenvector[i] = new double[prior_total_on];
	}
	double eigenvalue[prior_total_on];

	prior_alloc_subcovar(data, subcovar, eigenvector, eigenvalue, prior_subcovar_index);
	data.prior_eigenvector = eigenvector;
	data.prior_eigenvalue = eigenvalue;
	data.prior_subcovar = subcovar;
	data.prior_subcovar_index = prior_subcovar_index;


	data.fit_galaxy = false;
	shutoff_global(pars, data);

//	shutoff_sn_pars(pars, data);

	set_star_coef(a, data, 10, -0.4);
//  Sets Cardelli et al. (1989) reddening law with R_V = 3.1. Should not be changed, because other global parameters depend on it.
	set_cardelli(a, data, 3.1);

	//return the best-fit parameters with cmpfit errorbars
	double uncert[n_coef], uncert_rescaled[n_coef];
	double covar[(long)n_coef*n_coef], covar_dummy[n_coef*n_coef];
	for (long i=0;i<(n_coef*n_coef);i++) covar[i] = 0.0;
	for (int i=0;i<1000;i++) for(int j=0;j<1000;j++) global_covar[i][j] = 0.0;

	int nfree_tot=0;
	for (int i=0;i<n_coef;i++) {
		//printf("%i %f %i\n", i, a[i], pars[i].fixed);
		nfree_tot += 1-pars[i].fixed;
	}
	printf("nfree_tot: %i\n", nfree_tot);

	int nfitdata = 0; bool blocked;
	for (int i=0;i<data.n_data;i++) {
		int star_start = data.n_per_filter*data.n_phot_bands+data.n_velo_lines+data.n_gal+data.n_per_star*data.name[i];
		if (data.fit_star[data.name[i]] == 0) {continue;}
		if ((data.hjd[i]-a[star_start+0]) < 0) {continue;}
		if ((data.hjd[i]-a[star_start+0]) > 500) {continue;}
		if (data.dset[i] == 0) {
			if (data.flt[i] != 0) {continue;}
			if (data.fit_vel == 0) {continue;}
			if (data.fit_vel_star[data.name[i]] == 0) {continue;}
		}
		nfitdata++;
}
	printf("nfitdata: %i\n", nfitdata);


	// identify parameters inconsistent with constraints
	    for (int i=0; i<n_coef; i++) {
	      if ( (pars[i].limited[0] && (a[i] < pars[i].limits[0])) ||  (pars[i].limited[1] && (a[i] > pars[i].limits[1])) ) {
			int nstar = (i-data.n_per_filter*data.n_phot_bands - data.n_velo_lines - data.n_gal)/data.n_per_star;
			printf("inconsistency: %i %i %f  %f lower: %f upper: %f\n", i, nstar, a[i], data.min_time[nstar], pars[i].limits[0], pars[i].limits[1]);
			return 0;
      }}

	double chi2;
	chi2 = perform_fit(data, pars, a, uncert, covar, n_coef, 0);

	printf("Calculating bolometric light curve (lbol_err.dat)...\n");
	calculate_lbol(data, a);
	double mni, mni_err, l50, l50_err;
	//calculate_mni(data, a, mni, mni_err, l50, l50_err);  // this function does not compile

	write_data(a,data);
	write_model(a,data);

	//printf("calculating ln85...\n"); //values not rescaled
	//calculate_ln85(data,a,0);
	//calculate_popov93(data,a,0);

	// renormalize data errorbars to get chi^2/DOF = 1
	for (int i=0;i<data.n_data;i++) {
		data.err[i] *= sqrt(chi2/(nfitdata-nfree_tot));}

        printf("calculating ln85 & popov93...\n"); //values rescaled
	calculate_ln85(data,a,0);
	calculate_popov93(data,a,0);
	
	printf("chi2 with original errorbars: %f\n", chi2);
	for (int i=0;i<n_coef;i++) a_rescale[i] = a[i];
	chi2 = perform_fit(data, pars, a_rescale, uncert_rescaled, covar_dummy, n_coef, 0);
	printf("chi2 with rescaled errorbars (chi2/DOF=1): %f\n", chi2);

	printf("calculating mni_dependence...\n");
	mni_dependence(data,a,0);

        // obtain mni and lbol (this is not the best way of doing it)
        string str_logmni;	// log(Mni)
	string str_logmni_var;  // var(log_Mni)
        string str_l50;		// log(Lbol) @ 50 days
	string str_l50_var;	//var(log_Lbol)
	string skip;
	ifstream myfile( "popov93_err_rescale.dat" );
    	std::getline(myfile, skip, ' '); //skip the first value
    	std::getline(myfile, skip, ' '); //skip the second value.... etc  (check ln85.cpp for more info about "popov93_err_rescale.dat")
    	std::getline(myfile, skip, ' '); 
	std::getline(myfile, str_logmni, ' ');
	std::getline(myfile, str_l50, ' ');
    	std::getline(myfile, skip, ' ');
    	std::getline(myfile, skip, ' ');
    	std::getline(myfile, skip, ' ');
	std::getline(myfile, str_logmni_var,' ');
    	std::getline(myfile, skip, ' ');
    	std::getline(myfile, skip, ' ');
    	std::getline(myfile, skip, ' ');
    	std::getline(myfile, skip, ' ');
    	std::getline(myfile, skip, ' ');
	std::getline(myfile, str_l50_var,' ');

	mni = pow(10, ::atof(str_logmni.c_str()));	
	l50 = ::atof(str_l50.c_str());	
 	double s1 = ::atof(str_logmni_var.c_str());
	mni_err = mni*(sqrt(s1)*log(10));	// error propagation
	double s2 = ::atof(str_l50_var.c_str());
	l50_err = sqrt(s2);
  	//////////

	printf("******************\n");
	printf("Parameter Value Uncertainty Uncertainty_chi2/dof=1\n");
	printf("------------------\n");
        
	printf("Explosion time: %6.3f %5.3f %5.3f\n", a[star_start+0]+data.min_time[0], uncert[star_start+0], uncert_rescaled[star_start+0]);
	printf("Plateau duration: %6.3f %5.3f %5.3f\n", a[star_start+1], uncert[star_start+1], uncert_rescaled[star_start+1]);
	printf("Transition width:  %6.3f %5.3f %5.3f\n", a[star_start+2], uncert[star_start+2], uncert_rescaled[star_start+2]);
	printf("omega_0 [km/s]: %6.3f %5.1f %5.1f\n", a[star_start+3], uncert[star_start+3], uncert_rescaled[star_start+3]);
	printf("omega_1:        %7.4f %6.4f %6.4f\n", a[star_start+4], uncert[star_start+4], uncert_rescaled[star_start+4]);
	printf("omega_2 [km/s]: %6.1f %5.1f %5.1f\n", a[star_start+5], uncert[star_start+5], uncert_rescaled[star_start+5]);
	printf("gamma_0: 	%8.5f %7.5f %7.5f\n", a[star_start+6], uncert[star_start+6], uncert_rescaled[star_start+6]);
	printf("gamma_1: %6.3f %5.3f %5.3f\n", a[star_start+7], uncert[star_start+7], uncert_rescaled[star_start+7]);
	printf("alpha_0: %8.5f %7.5f %7.5f\n", a[star_start+8], uncert[star_start+8], uncert_rescaled[star_start+8]);
	printf("alpha_1: %6.4f %5.4f %5.4f\n", a[star_start+9], uncert[star_start+9], uncert_rescaled[star_start+9]);
	printf("Distance modulus: %6.3f %5.3f %5.3f\n", a[star_start+11]+30, uncert[star_start+11], uncert_rescaled[star_start+11]);
	printf("E(B-V): %6.4f %5.4f %5.4f\n", a[star_start+12], uncert[star_start+12], uncert_rescaled[star_start+12]);

	printf("------------------\n");
	printf("Derived parameters\n");
	printf("------------------\n");

	printf("Ni mass [Msun]: %5.3f +/- %5.3f\n", mni, mni_err);
	if ((data.max_time[0]-data.min_time[0]-a[0]) < 200) printf("Warning: Ni mass is an extrapolation, there are no datapoints more than 200 days after explosion!\n");
	printf("Lbol@50 days [Lsun]: %10.3le  +/- %10.3le\n", l50, l50_err);
	if (a[0] < -50) printf("Warning: Lbol@50 days is an extrapolation, there are no datapoints less than 50 days after explosion!\n");
	printf("log E_exp/10^50 ergs (LN85): %f\n", ln85_e(data, a_rescale, 0, 2, 0));
	printf("log M_ej/Msun (LN85): %f\n", ln85_m(data, a_rescale, 0, 2, 0));
	printf("log R_prog/Rsun (LN85): %f\n", ln85_r(data, a_rescale, 0, 2, 0));
	printf("log E_exp/10^50 ergs (POPOV93): %f\n", popov93_e(data, a_rescale, 0, 2, 0));
	printf("log M_ej/Msun (POPOV93): %f\n", popov93_m(data, a_rescale, 0, 2, 0));
	printf("log R_prog/Rsun (POPOV93): %f\n", popov93_r(data, a_rescale, 0, 2, 0));

        
        // here I modify the output files names to include the SN name given at the beginning
        string phot = sn_name + "/" + sn_name + "_phot.dat";
	string model_phot = sn_name + "/" + sn_name + "_model_phot.dat";
	string velo = sn_name + "/" + sn_name + "_velo.dat";
	string model_velo = sn_name + "/" + sn_name + "_model_velo.dat";
	string lbol_err = sn_name + "/" + sn_name + "_lbol_err.dat";
	string ln85_err_norescale = sn_name + "/" + sn_name + "_ln85_err_rescale.dat";
	string popov93_err_norescale = sn_name + "/" + sn_name + "_popov93_err_rescale.dat";
	
        std::ifstream  src_phot("phot.dat", std::ios::binary);
     	std::ofstream  dst_phot(phot.c_str(),   std::ios::binary);
     	dst_phot << src_phot.rdbuf();
	src_phot.close();
	dst_phot.close();
	
	std::ifstream  src_model_phot("model_phot.dat", std::ios::binary);
     	std::ofstream  dst_model_phot(model_phot.c_str(),   std::ios::binary);
     	dst_model_phot << src_model_phot.rdbuf();
	src_model_phot.close();
	dst_model_phot.close();

	std::ifstream  src_velo("velo.dat", std::ios::binary);
     	std::ofstream  dst_velo(velo.c_str(),   std::ios::binary);
     	dst_velo << src_velo.rdbuf();
	src_velo.close();
	dst_velo.close();

	std::ifstream  src_model_velo("model_velo.dat", std::ios::binary);
     	std::ofstream  dst_model_velo(model_velo.c_str(),   std::ios::binary);
     	dst_model_velo << src_model_velo.rdbuf();
	src_model_velo.close();
	dst_model_velo.close();

	std::ifstream  src_lbol_err("lbol_err.dat", std::ios::binary);
     	std::ofstream  dst_lbol_err(lbol_err.c_str(),   std::ios::binary);
     	dst_lbol_err << src_lbol_err.rdbuf();
	src_lbol_err.close();
	dst_lbol_err.close();

	std::ifstream  src_ln85("ln85_err_rescale.dat", std::ios::binary);
     	std::ofstream  dst_ln85(ln85_err_norescale.c_str(),   std::ios::binary);
     	dst_ln85 << src_ln85.rdbuf();
	src_ln85.close();
	dst_ln85.close();

	std::ifstream  src_popov93("popov93_err_rescale.dat", std::ios::binary);
     	std::ofstream  dst_popov93(popov93_err_norescale.c_str(),   std::ios::binary);
     	dst_popov93 << src_popov93.rdbuf();
	src_popov93.close();
	dst_popov93.close();
        
}
