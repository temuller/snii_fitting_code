/* deriv/deriv.c
 *
 * Copyright (C) 2004, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "var_struct.h"

#define GSL_DBL_EPSILON        2.2204460492503131e-16



static void central_deriv (double f(vars_struct, double *, const int, const int, const double), double *x, const int n, const int k, double h, vars_struct data, const int star_no, const int phot_band, const double time, double *result, double *abserr_round, double *abserr_trunc) {
  /* Compute the derivative using the 5-point rule (x-h, x-h/2, x,
     x+h/2, x+h). Note that the central point is not used.

     Compute the error using the difference between the 5-point and
     the 3-point rule (x-h,x,x+h). Again the central point is not
     used. */

  double xminus[n], xplus[n];
	for (int i=0;i<n;i++) {
		xminus[i] = x[i];
		xplus[i]  = x[i];
	}
	xminus[k] -= h;
	xplus[k] += h;

  double fm1 = f(data, xminus, star_no, phot_band, time);
  double fp1 = f(data, xplus, star_no, phot_band, time);

	for (int i=0;i<n;i++) {
		xminus[i] = x[i];
		xplus[i]  = x[i];
	}
	xminus[k] -= h/2;
	xplus[k] += h/2;

  double fmh = f(data, xminus, star_no, phot_band, time);
  double fph = f(data, xplus, star_no, phot_band, time);

  double r3 = 0.5 * (fp1 - fm1);
  double r5 = (4.0 / 3.0) * (fph - fmh) - (1.0 / 3.0) * r3;

  double e3 = (fabs (fp1) + fabs (fm1)) * GSL_DBL_EPSILON;
  double e5 = 2.0 * (fabs (fph) + fabs (fmh)) * GSL_DBL_EPSILON + e3;

  /* The next term is due to finite precision in x+h = O (eps * x) */

  double dy = fmax (fabs (r3 / h), fabs (r5 / h)) *(fabs (x[k]) / h) * GSL_DBL_EPSILON;

  /* The truncation error in the r5 approximation itself is O(h^4).
     However, for safety, we estimate the error from r5-r3, which is
     O(h^2).  By scaling h we will minimise this estimated error, not
     the actual truncation error in r5. */

  *result = r5 / h;
  *abserr_trunc = fabs ((r5 - r3) / h); /* Estimated truncation error O(h^2) */
  *abserr_round = fabs (e5 / h) + dy;   /* Rounding error (cancellations) */
}

double gsl_deriv_central (double f(vars_struct, double *, const int, const int, const double), double *x,const int n, const int k, double h, vars_struct data, const int star_no, const int phot_band, const double time) {
  double r_0, round, trunc, error;
  central_deriv (f, x, n, k, h,data, star_no, phot_band, time, &r_0, &round, &trunc);
  error = round + trunc;

  if (round < trunc && (round > 0 && trunc > 0))
    {
      double r_opt, round_opt, trunc_opt, error_opt;

      /* Compute an optimised stepsize to minimize the total error,
         using the scaling of the truncation error (O(h^2)) and
         rounding error (O(1/h)). */

      double h_opt = h * pow (round / (2.0 * trunc), 1.0 / 3.0);
      central_deriv (f, x,n,k, h_opt,data, star_no, phot_band, time, &r_opt, &round_opt, &trunc_opt);
      error_opt = round_opt + trunc_opt;

      /* Check that the new error is smaller, and that the new derivative
         is consistent with the error bounds of the original estimate. */

      if (error_opt < error && fabs (r_opt - r_0) < 4.0 * error)
        {
          r_0 = r_opt;
          error = error_opt;
        }
    }

//  *result = r_0;
//  *abserr = error;

  return r_0;
}