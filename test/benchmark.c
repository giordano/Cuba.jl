/* benchmark.c --- Benchmark Cuba Library in C
 *
 * Copyright (C) 2016  Mosè Giordano
 *
 * Maintainer: Mosè Giordano <mose AT gnu DOT org>
 * Keywords: numeric integration, benchmark
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Commentary:
 *
 * Based on file demo/demo-c.c of Cuba Library, by Thomas Hahn, released with
 * the same license.
 *
 * Code: */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cuba.h"

#define NDIM 3
#define NCOMP 11
#define USERDATA NULL
#define NVEC 1
#define EPSREL 1e-8
#define EPSABS 1e-8
#define VERBOSE 0
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 1000000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.0

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

static int Integrand(const int *ndim, const double xx[],
		     const int *ncomp, double ff[], void *userdata) {
#define x xx[0]
#define y xx[1]
#define z xx[2]
#define Sq(x) x*x
#define rsq (Sq(x) + Sq(y) + Sq(z))
  ff[0]  = sin(x)*cos(y)*exp(z);
  ff[1]  = 1./((x + y)*(x + y) + .003)*cos(y)*exp(z);
  ff[2]  = 1/(3.75 - cos(M_PI*x) - cos(M_PI*y) - cos(M_PI*z));
  ff[3]  = fabs(rsq - .125);
  ff[4]  = exp(-rsq);
  ff[5]  = 1./(1. - x*y*z + 1e-10);
  ff[6]  = sqrt(fabs(x - y - z));
  ff[7]  = exp(-x*y*z);
  ff[8]  = Sq(x)/(cos(x + y + z + 1) + 5);
  ff[9]  = (x > .5) ? 1/sqrt(x*y*z + 1e-5) : sqrt(x*y*z);
  ff[10] = (rsq < 1) ? 1 : 0;
  return 0;
}

int main() {
  int comp, nregions, fail;
  long long int neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
  float start, end;
  start = (float)clock();
  llVegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
	  EPSREL, EPSABS, VERBOSE, SEED,
	  MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
	  GRIDNO, STATEFILE, SPIN,
	  &neval, &fail, integral, error, prob);
  end = (float)clock();
  printf("%10.6f seconds (Vegas)\n", (end - start)/CLOCKS_PER_SEC);
  start = (float)clock();
  llSuave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
	  EPSREL, EPSABS, VERBOSE | LAST, SEED,
	  MINEVAL, MAXEVAL, NNEW, NMIN, 25.0,
	  STATEFILE, SPIN,
	  &nregions, &neval, &fail, integral, error, prob);
  end = (float)clock();
  printf("%10.6f seconds (Suave)\n", (end - start)/CLOCKS_PER_SEC);
  start = (float)clock();
  llDivonne(NDIM, NCOMP, Integrand, USERDATA, NVEC,
	    EPSREL, EPSABS, VERBOSE, SEED,
	    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
	    BORDER, MAXCHISQ, MINDEVIATION,
	    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
	    STATEFILE, SPIN,
	    &nregions, &neval, &fail, integral, error, prob);
  end = (float)clock();
  printf("%10.6f seconds (Divonne)\n", (end - start)/CLOCKS_PER_SEC);
  start = (float)clock();
  llCuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
	  EPSREL, EPSABS, VERBOSE | LAST,
	  MINEVAL, MAXEVAL, KEY,
	  STATEFILE, SPIN,
	  &nregions, &neval, &fail, integral, error, prob);
  end = (float)clock();
  printf("%10.6f seconds (Cuhre)\n", (end - start)/CLOCKS_PER_SEC);
  return 0;
}
