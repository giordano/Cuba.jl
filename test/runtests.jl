### runtests.jl --- Test suite for Cuba.jl

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

### Code:

using Cuba
using Base.Test

f1(x,y,z) = sin(x)*cos(y)*exp(z)
f2(x,y,z) = exp(-(x*x + y*y + z*z))
f3(x,y,z) = 1/(1 - x*y*z)

# This is the integrand function.  This is just like one would do in C.
function func(ndim::Cint,
              xx::Ptr{Cdouble},
              ncomp::Cint,
              ff::Ptr{Cdouble},
              userdata::Ptr{Void}=USERDATA_DEF
              )
    x = pointer_to_array(xx, (ndim,))
    f = pointer_to_array(ff, (ncomp,))

    # Define components of the integrand function.
    f[1] = f1(x[1], x[2], x[3])
    f[2] = f2(x[1], x[2], x[3])
    f[3] = f3(x[1], x[2], x[3])

    xx = pointer_from_objref(x)
    ff = pointer_from_objref(f)
    return Cint(0)::Cint
end

# Accelerators and cores settings
Cuba.accel(0,1000)
Cuba.cores(0,1000)

# Test results and make sure the estimation of error is exact.
let
    local result, res1, res2, res3
    res1 = (e-1)*(1-cos(1))*sin(1)
    res2 = (sqrt(pi)*erf(1)/2)^3
    res3 = zeta(3)
    # Vegas
    result = Vegas(func, ndim=3, ncomp=3, epsabs=1e-4, epsrel=1e-8)
    @test_approx_eq_eps result[1][1]  res1  result[2][1]
    @test_approx_eq_eps result[1][2]  res2  result[2][2]
    @test_approx_eq_eps result[1][3]  res3  result[2][3]
    # Suave
    result = Suave(func, ndim=3, ncomp=3, epsabs=1e-3, epsrel=1e-8)
    @test_approx_eq_eps result[1][1]  res1  result[2][1]
    @test_approx_eq_eps result[1][2]  res2  result[2][2]
    @test_approx_eq_eps result[1][3]  res3  result[2][3]
    # Divonne
    result = Divonne(func, ndim=3, ncomp=3, epsabs=1e-4, epsrel=1e-8)
    @test_approx_eq_eps result[1][1]  res1  result[2][1]
    @test_approx_eq_eps result[1][2]  res2  result[2][2]
    # @test_approx_eq_eps result[1][3]  res3  result[2][3] # <== This integral diverges!
    # Cuhre
    result = Cuhre(func, ndim=3, ncomp=3, epsabs=1e-8, epsrel=1e-8)
    @test_approx_eq_eps result[1][1]  res1  result[2][1]
    @test_approx_eq_eps result[1][2]  res2  result[2][2]
    @test_approx_eq_eps result[1][3]  res3  result[2][3]
end

# Test other function
Cuba.init(C_NULL, C_NULL)
Cuba.exit(C_NULL, C_NULL)
