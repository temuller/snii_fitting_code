#include <stdlib.h>
#include <math.h>


double reddening_coef(double lambda, double rv) {
	double x = 1.0/lambda;

if (x <= 1.1) {
	double a1, b1;
	a1 =  0.574*pow(x, 1.61);
	b1 = -0.527*pow(x, 1.61);
	return a1*rv + b1;
}

if ((x > 1.1) && (x <= 3.3)) {
	double y = x-1.82;
	double a2 = 1.0+0.17699*y-0.50447*y*y-0.02427*pow(y,3)+0.72085*pow(y,4) + 0.01979*pow(y,5)-0.77530*pow(y,6)+0.32999*pow(y,7);
	double b2 = 1.41338*y+2.28305*y*y+1.07233*pow(y,3)-5.38434*pow(y,4)-0.62251*pow(y,5)+5.30260*pow(y,6)-2.09002*pow(y,7);
	return a2*rv + b2;
}


if ( (x > 3.3) && (x <= 8)) {
	double fa = 0, fb = 0;

	if (x >= 5.9) {
		fa = -0.04473*pow(x-5.9, 2)-0.009779*pow(x-5.9,3);
		fb = 0.2130*pow(x-5.9,2) + 0.1207*pow(x-5.9,3);
	}
	double a3 = 1.752 - 0.316*x - 0.104/(pow(x-4.67,2) + 0.341) + fa;
	double b3 = -3.090 + 1.825*x + 1.203/(pow(x-4.62,2) + 0.263) + fb;
	return a3*rv + b3;
}

if ( (x > 8) && (x <= 10) ) {
	double a4 = -1.073 - 0.628*(x-8) + 0.137*pow(x-8, 2) - 0.070*pow(x-8,3);
	double b4 = 13.670 + 4.257*(x-8) - 0.420*pow(x-8, 2) + 0.374*pow(x-8,3);
	return a4*rv + b4;
}
return 0;

}
