c     benchmark.f --- Benchmark Cuba Library in Fortran
c
c     Copyright (C) 2016  Mosè Giordano
c
c     Maintainer: Mosè Giordano <mose AT gnu DOT org>
c     Keywords: numeric integration, benchmark
c
c     This program is free software: you can redistribute it and/or
c     modify it under the terms of the GNU Lesser General Public License
c     as published by the Free Software Foundation, either version 3 of
c     the License, or (at your option) any later version.
c
c     This program is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     Lesser General Public License for more details.
c
c     You should have received a copy of the GNU Lesser General Public
c     License along with this program.  If not, see
c     <http://www.gnu.org/licenses />.
c
c     Commentary:
c
c     Based on file demo/demo-fortran.F of Cuba Library, by Thomas Hahn,
c     released with the same license.
c
c     Code:

      program CubaTest
      implicit none

      integer ndim, ncomp, last, seed
      integer(8) nvec, mineval, maxeval
      double precision epsrel, epsabs, userdata
      parameter (ndim = 3)
      parameter (ncomp = 11)
      parameter (userdata = 0)
      parameter (nvec = 1)
      parameter (epsrel = 1D-8)
      parameter (epsabs = 1D-8)
      parameter (last = 4)
      parameter (seed = 0)
      parameter (mineval = 0)
      parameter (maxeval = 1000000)

      integer(8) nstart, nincrease, nbatch, gridno
      integer*8 spin
      character*(*) statefile
      parameter (nstart = 1000)
      parameter (nincrease = 500)
      parameter (nbatch = 1000)
      parameter (gridno = 0)
      parameter (statefile = "")
      parameter (spin = -1)

      integer(8) nnew, nmin
      double precision flatness
      parameter (nnew = 1000)
      parameter (nmin = 2)
      parameter (flatness = 25D0)

      integer key1, key2, key3, maxpass
      double precision border, maxchisq, mindeviation
      integer ldxgiven
      integer(8) ngiven, nextra
      parameter (key1 = 47)
      parameter (key2 = 1)
      parameter (key3 = 1)
      parameter (maxpass = 5)
      parameter (border = 0D0)
      parameter (maxchisq = 10D0)
      parameter (mindeviation = .25D0)
      parameter (ngiven = 0)
      parameter (ldxgiven = ndim)
      parameter (nextra = 0)

      integer key
      parameter (key = 0)

      external integrand

      double precision integral(ncomp), error(ncomp), prob(ncomp)
      integer verbose, nregions, fail
      integer(8) neval

      character*16 env
      parameter (verbose = 0)

      real start, finish

      integer c

      call cpu_time(start)
      call llvegas(ndim, ncomp, integrand, userdata, nvec,
     &    epsrel, epsabs, verbose, seed,
     &    mineval, maxeval, nstart, nincrease, nbatch,
     &    gridno, statefile, spin,
     &    neval, fail, integral, error, prob)
      call cpu_time(finish)
      print '(f10.6, a)', finish-start, " seconds (Vegas)"

      call cpu_time(start)
      call llsuave(ndim, ncomp, integrand, userdata, nvec,
     &    epsrel, epsabs, verbose + last, seed,
     &    mineval, maxeval, nnew, nmin, flatness,
     &    statefile, spin,
     &    nregions, neval, fail, integral, error, prob)
      call cpu_time(finish)
      print '(f10.6, a)', finish-start, " seconds (Suave)"

      call cpu_time(start)
      call lldivonne(ndim, ncomp, integrand, userdata, nvec,
     &    epsrel, epsabs, verbose, seed,
     &    mineval, maxeval, key1, key2, key3, maxpass,
     &    border, maxchisq, mindeviation,
     &    ngiven, ldxgiven, 0, nextra, 0,
     &    statefile, spin,
     &    nregions, neval, fail, integral, error, prob)
      call cpu_time(finish)
      print '(f10.6, a)', finish-start, " seconds (Divonne)"

      call cpu_time(start)
      call llcuhre(ndim, ncomp, integrand, userdata, nvec,
     &    epsrel, epsabs, verbose + last,
     &    mineval, maxeval, key,
     &    statefile, spin,
     &    nregions, neval, fail, integral, error, prob)
      call cpu_time(finish)
      print '(f10.6, a)', finish-start, " seconds (Cuhre)"

      end

      integer function integrand(ndim, xx, ncomp, ff)
      implicit none
      integer ndim, ncomp
      double precision xx(*), ff(*)

#define x xx(1)
#define y xx(2)
#define z xx(3)

      double precision pi, rsq
      parameter (pi = 3.14159265358979323846D0)

      rsq = x**2 + y**2 + z**2

      ff(1) = sin(x)*cos(y)*exp(z)
      ff(2) = 1/((x + y)**2 + .003D0)*cos(y)*exp(z)
      ff(3) = 1/(3.75D0 - cos(pi*x) - cos(pi*y) - cos(pi*z))
      ff(4) = abs(rsq - .125D0)
      ff(5) = exp(-rsq)
      ff(6) = 1/(1 - x*y*z + 1D-10)
      ff(7) = sqrt(abs(x - y - z))
      ff(8) = exp(-x*y*z)
      ff(9) = x**2/(cos(x + y + z + 1) + 5)
      if( x .gt. .5D0 ) then
        ff(10) = 1/sqrt(x*y*z + 1D-5)
      else
        ff(10) = sqrt(x*y*z)
      endif
      if( rsq .lt. 1 ) then
        ff(11) = 1
      else
        ff(11) = 0
      endif
      integrand = 0
      end
