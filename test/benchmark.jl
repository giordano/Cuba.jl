### benchmark.jl --- Benchmark Cuba.jl and Cuba C Library

# Copyright (C) 2016  Mosè Giordano

# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: numeric integration

# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

### Commentary:

# Load this file in order to compare performance of Cuba.jl and original Cuba
# Library in C.

### Code:

using Cuba

ndim=3
ncomp=11
epsabs=1e-8
epsrel=1e-8

Sq(x)      = x*x
rsq(x,y,z) = Sq(x) + Sq(y) + Sq(z)
t1(x,y,z)  = sin(x)*cos(y)*exp(z)
t2(x,y,z)  = 1./((x + y)*(x + y) + .003)*cos(y)*exp(z)
t3(x,y,z)  = 1/(3.75 - cos(pi*x) - cos(pi*y) - cos(pi*z))
t4(x,y,z)  = abs(rsq(x,y,z) - .125)
t5(x,y,z)  = exp(-rsq(x,y,z))
t6(x,y,z)  = 1./(1. - x*y*z + 1e-10)
t7(x,y,z)  = sqrt(abs(x - y - z))
t8(x,y,z)  = exp(-x*y*z)
t9(x,y,z)  = Sq(x)/(cos(x + y + z + 1) + 5)
t10(x,y,z) = (x > .5) ? 1/sqrt(x*y*z + 1e-5) : sqrt(x*y*z)
t11(x,y,z) = (rsq(x,y,z) < 1) ? 1 : 0
function test(ndim::Cint, xx::Ptr{Cdouble}, ncomp::Cint,
              ff::Ptr{Cdouble}, u::Ptr{Void}=C_NULL)
    x = pointer_to_array(xx, (ndim,))
    f = pointer_to_array(ff, (ncomp,))
    f[1]  = t1( x[1], x[2], x[3])
    f[2]  = t2( x[1], x[2], x[3])
    f[3]  = t3( x[1], x[2], x[3])
    f[4]  = t4( x[1], x[2], x[3])
    f[5]  = t5( x[1], x[2], x[3])
    f[6]  = t6( x[1], x[2], x[3])
    f[7]  = t7( x[1], x[2], x[3])
    f[8]  = t8( x[1], x[2], x[3])
    f[9]  = t9( x[1], x[2], x[3])
    f[10] = t10(x[1], x[2], x[3])
    f[11] = t11(x[1], x[2], x[3])
    ff = pointer_from_objref(f)
    return Cint(0)::Cint
end

info("Ignore these times...")
@time Vegas(test, ndim=ndim, ncomp=ncomp, epsabs=epsabs, epsrel=epsrel);
@time Suave(test, ndim=ndim, ncomp=ncomp, epsabs=epsabs, epsrel=epsrel);
@time Divonne(test, ndim=ndim, ncomp=ncomp, epsabs=epsabs, epsrel=epsrel);
@time Cuhre(test, ndim=ndim, ncomp=ncomp, epsabs=epsabs, epsrel=epsrel);

info("Performance of Cuba.jl:")
@time Vegas(test, ndim=ndim, ncomp=ncomp, epsabs=epsabs, epsrel=epsrel);
@time Suave(test, ndim=ndim, ncomp=ncomp, epsabs=epsabs, epsrel=epsrel);
@time Divonne(test, ndim=ndim, ncomp=ncomp, epsabs=epsabs, epsrel=epsrel);
@time Cuhre(test, ndim=ndim, ncomp=ncomp, epsabs=epsabs, epsrel=epsrel);

cd(dirname(@__FILE__))
cp("../deps/cuba-shared-object/libcuba.a", "libcuba.a", remove_destination=true)
cp("../deps/cuba-shared-object/cuba.h", "cuba.h", remove_destination=true)
run(`gcc -o benchmark benchmark.c libcuba.a -lm`)
info("Performance of Cuba C Library:")
run(`./benchmark`)
