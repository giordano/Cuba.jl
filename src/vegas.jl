### vegas.jl --- Integrate with Vegas method

# Copyright (C) 2017  Mosè Giordano

# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: numerical integration

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

immutable Vegas{T, I} <: Integral{T, I}
    func::T
    ndim::Int
    ncomp::Int
    nvec::I
    reltol::Cdouble
    abstol::Cdouble
    flags::Int
    seed::Int
    minevals::I
    maxevals::I
    nstart::I
    nincrease::I
    nbatch::I
    gridno::Int
    statefile::String
    spin::Ptr{Void}
end

for (CubaInt, prefix) in ((Int32, ""), (Int64, "ll"))
    @eval @inline function dointegrate!{T}(x::Vegas{T, $CubaInt}, integrand, nregions,
                                           neval, fail, integral, error, prob)
        ccall(($(prefix * "Vegas"), libcuba), Cdouble,
              (Cint, # ndim
               Cint, # ncomp
               Ptr{Void}, # integrand
               Any, # userdata
               $CubaInt, # nvec
               Cdouble, # reltol
               Cdouble, # abstol
               Cint, # flags
               Cint, # seed
               $CubaInt, # minevals
               $CubaInt, # maxevals
               $CubaInt, # nstart
               $CubaInt, # nincrease
               $CubaInt, # nbatch
               Cint, # gridno
               Ptr{Cchar}, # statefile
               Ptr{Void}, # spin
               Ptr{$CubaInt}, # neval
               Ptr{Cint}, # fail
               Ptr{Cdouble}, # integral
               Ptr{Cdouble}, # error
               Ptr{Cdouble}),# prob
              # Input
              x.ndim, x.ncomp, integrand, x.func, x.nvec,
              x.reltol, x.abstol, x.flags, x.seed, x.minevals, x.maxevals,
              x.nstart, x.nincrease, x.nbatch, x.gridno, x.statefile, x.spin,
              # Output
              neval, fail, integral, error, prob)
    end

    func = Symbol(prefix, "vegas")
    @eval function $func{T}(integrand::T, ndim::Integer=1, ncomp::Integer=1;
                            nvec::Integer=NVEC, reltol::Real=RELTOL, abstol::Real=ABSTOL,
                            flags::Integer=FLAGS, seed::Integer=SEED,
                            minevals::Real=MINEVALS, maxevals::Real=MAXEVALS,
                            nstart::Integer=NSTART, nincrease::Integer=NINCREASE,
                            nbatch::Integer=NBATCH, gridno::Integer=GRIDNO,
                            statefile::AbstractString=STATEFILE, spin::Ptr{Void}=SPIN)
        return dointegrate(Vegas(integrand, ndim, ncomp, $CubaInt(nvec), Cdouble(reltol),
                                 Cdouble(abstol), flags, seed, trunc($CubaInt, minevals),
                                 trunc($CubaInt, maxevals), $CubaInt(nstart),
                                 $CubaInt(nincrease), $CubaInt(nbatch), gridno,
                                 String(statefile), spin))
    end
end

"""
    vegas(integrand, ndim=1, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions
    llvegas(integrand, ndim=1, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions

Calculate integral of `integrand` over the unit hypercube in `ndim` dimensions
using Vegas algorithm.  `integrand` is a vectorial function with `ncomp`
components.  `ndim` and `ncomp` default to 1.

Accepted keywords:

* `nvec`
* `reltol`
* `abstol`
* `flags`
* `seed`
* `minevals`
* `maxevals`
* `nstart`
* `nincrease`
* `nbatch`
* `gridno`
* `statefile`
* `spin`
"""
vegas, llvegas
