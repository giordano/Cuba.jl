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
function func(ndim::Cint, xx::Ptr{Cdouble}, ncomp::Cint, ff::Ptr{Cdouble},
              userdata::Ptr{Void})
    # Take arrays from "xx" and "ff" pointers.
    x = pointer_to_array(xx, (ndim,))
    f = pointer_to_array(ff, (ncomp,))
    # Define components of the integrand function.
    f[1] = f1(x[1], x[2], x[3])
    f[2] = f2(x[1], x[2], x[3])
    f[3] = f3(x[1], x[2], x[3])
    # Store back the results to "ff"
    ff = pointer_from_objref(f)
    return Cint(0)::Cint
end

# Accelerators and cores settings
Cuba.accel(0,1000)
Cuba.cores(0,1000)

# Test results and make sure the estimation of error is exact.
epsabs = Dict(Vegas=>1e-4, Suave=>1e-3, Divonne=>1e-4, Cuhre=>1e-8)
answer = [(e-1)*(1-cos(1))*sin(1), (sqrt(pi)*erf(1)/2)^3, zeta(3)]
ncomp = 3
for alg in (:Vegas, :Suave, :Divonne, :Cuhre)
    info("Testing $(string(alg)) algorithm")
    result = @eval $alg($func, 3, $ncomp, epsabs=$epsabs[$alg],
                        epsrel=1e-8)
    for i = 1:ncomp
        println("Component $i: ", result[1][i], " ± ", result[2][i],
                isfinite(result[1][i]) ? "" : " (skipping test)")
        println("Should be:   ", answer[i])
        if isfinite(result[1][i])
            @test_approx_eq_eps result[1][i] answer[i] result[2][i]
        end
    end
end

# Test other function
Cuba.init(C_NULL, C_NULL)
Cuba.exit(C_NULL, C_NULL)
