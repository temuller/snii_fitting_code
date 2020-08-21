#include "math.h"

void define_filter_zeropoints(double *fltabs, double *lam, const int n_phot_bands) {
	for (int i=0;i<n_phot_bands;i++) fltabs[i] = 0.0;
	const double c =2.99792458e10;
	// now for UBVRIJHKLM from Johnson 1965
	fltabs[0] = 3.98e-9; // 4.35e-12*1e7*1e-4; //3.98e-9; // U, erg cm^-2 s^-1 A^-1
	fltabs[1] = 6.95e-9; //B7.20e-12*1e7*1e-4; //6.95e-9; //B
	fltabs[2] = 3.63e-9; //V3.92e-12*1e7*1e-4; //3.63e-9; //V
	fltabs[3] = 2.254e-9; //Rc1.76e-12*1e7*1e-4; //2.254e-9; //Rc
	fltabs[4] = 1.196e-9; //Ic 8.3e-13*1e7*1e-4;  //1.196e-9; //Ic
	fltabs[5] = 3.4e-13*1e7*1e-4;  //W cm^-2 mikro m^-1 Johnson JHK
	fltabs[6] = 1.26e-13*1e7*1e-4; //"
	fltabs[7] = 3.9e-14*1e7*1e-4; //"
	fltabs[8] = 3631*1e-23*c/(   (lam[8]*1e-4) * (lam[8]*1e4) )   ;//SDSS conversion from Jansky to Flambda
	fltabs[9] =  3631*1e-23*c/(   (lam[9]*1e-4) * (lam[9]*1e4) ) ;
	fltabs[10] = 3631*1e-23*c/(   (lam[10]*1e-4) * (lam[10]*1e4) );
	fltabs[11] = 3631*1e-23*c/(   (lam[11]*1e-4) * (lam[11]*1e4) );
	fltabs[12] = 3631*1e-23*c/(   (lam[12]*1e-4) * (lam[12]*1e4) );
	fltabs[13] = pow(10, 0.4*17.35)*6.0e-16; //Swift bands from Poole et al. (2008)
	fltabs[14] = pow(10, 0.4*16.82)*7.5e-16;
	fltabs[15] = pow(10, 0.4*17.49)*4.3e-16;
	fltabs[16] = pow(10, 0.4*18.34)*1.5e-16;
	fltabs[17] = pow(10, 0.4*19.11)*1.32e-16;
	fltabs[18] = pow(10, 0.4*17.89)*2.61e-16;
	fltabs[19] = 0;
	fltabs[20] = 0;
	fltabs[23] = 280.9*1e-23*(c/(lam[23]*1e-4))/(lam[23]*1e4); // Spitzer IRAC
	fltabs[24] = 179.7*1e-23*(c/(lam[24]*1e-4))/(lam[24]*1e4);
	fltabs[25] = 115.0*1e-23*(c/(lam[25]*1e-4))/(lam[25]*1e4);
	fltabs[26] = 64.9*1e-23*(c/(lam[26]*1e-4))/(lam[26]*1e4);
	fltabs[27] = 0;
}

