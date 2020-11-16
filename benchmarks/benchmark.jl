### benchmark.jl --- Benchmark Cuba.jl and Cuba C Library

# Copyright (C) 2016-2019  Mosè Giordano

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

using Cuba, Printf

const ndim=3
const ncomp=11
const atol=1e-8
const rtol=1e-8

rsq(x,y,z) = abs2(x) + abs2(y) + abs2(z)
t1(x,y,z)  = sin(x)*cos(y)*exp(z)
t2(x,y,z)  = 1.0/((x + y)*(x + y) + 0.003)*cos(y)*exp(z)
t3(x,y,z)  = 1.0/(3.75 - cos(pi*x) - cos(pi*y) - cos(pi*z))
t4(x,y,z)  = abs(rsq(x,y,z) - 0.125)
t5(x,y,z)  = exp(-rsq(x,y,z))
t6(x,y,z)  = 1.0/(1.0 - x*y*z + 1e-10)
t7(x,y,z)  = sqrt(abs(x - y - z))
t8(x,y,z)  = exp(-x*y*z)
t9(x,y,z)  = abs2(x)/(cos(x + y + z + 1.0) + 5.0)
t10(x,y,z) = (x > 0.5) ? 1.0/sqrt(x*y*z + 1e-5) : sqrt(x*y*z)
t11(x,y,z) = (rsq(x,y,z) < 1.0) ? 1.0 : 0.0
function test(x::Vector{Float64}, f::Vector{Float64})
    @inbounds f[1]  = t1( x[1], x[2], x[3])
    @inbounds f[2]  = t2( x[1], x[2], x[3])
    @inbounds f[3]  = t3( x[1], x[2], x[3])
    @inbounds f[4]  = t4( x[1], x[2], x[3])
    @inbounds f[5]  = t5( x[1], x[2], x[3])
    @inbounds f[6]  = t6( x[1], x[2], x[3])
    @inbounds f[7]  = t7( x[1], x[2], x[3])
    @inbounds f[8]  = t8( x[1], x[2], x[3])
    @inbounds f[9]  = t9( x[1], x[2], x[3])
    @inbounds f[10] = t10(x[1], x[2], x[3])
    @inbounds f[11] = t11(x[1], x[2], x[3])
end

@info "Performance of Cuba.jl:"
for alg in (vegas, suave, divonne, cuhre)
    # Run the integrator a first time to compile the function.
    alg(test, ndim, ncomp, atol=atol,
         rtol=rtol);
    start_time = time_ns()
    alg(test, ndim, ncomp, atol=atol,
        rtol=rtol);
    end_time = time_ns()
    println(@sprintf("%10.6f", Int(end_time - start_time)/1e9),
            " seconds (", uppercasefirst(string(nameof(alg))), ")")
end

cd(@__DIR__) do
    if mtime("benchmark.c") > mtime("benchmark-c")
        run(`gcc -O3 -I $(Cuba.Cuba_jll.artifact_dir)/include -o benchmark-c benchmark.c $(Cuba.Cuba_jll.libcuba_path) -lm`)
    end
    @info "Performance of Cuba Library in C:"
    withenv(Cuba.Cuba_jll.JLLWrappers.LIBPATH_env => Cuba.Cuba_jll.LIBPATH[]) do
        run(`./benchmark-c`)
    end

    if success(`which gfortran`)
        if mtime("benchmark.f") > mtime("benchmark-fortran")
            run(`gfortran -O3 -fcheck=no-bounds -cpp -o benchmark-fortran benchmark.f $(Cuba.Cuba_jll.libcuba_path) -lm`)
        end
        @info "Performance of Cuba Library in Fortran:"
        withenv(Cuba.Cuba_jll.JLLWrappers.LIBPATH_env => Cuba.Cuba_jll.LIBPATH[]) do
            run(`./benchmark-fortran`)
        end
    end
end
