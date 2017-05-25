### cuhre.jl --- Integrate with Cuhre method

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

immutable Cuhre{T, I} <: Integral{T, I}
    func::T
    ndim::Int
    ncomp::Int
    nvec::I
    reltol::Cdouble
    abstol::Cdouble
    flags::Int
    minevals::I
    maxevals::I
    key::Int
    statefile::String
    spin::Ptr{Void}
end

for (CubaInt, prefix) in ((Int32, ""), (Int64, "ll"))
    @eval @inline function dointegrate!{T}(x::Cuhre{T, $CubaInt}, integrand, nregions,
                                           neval, fail, integral, error, prob)
        ccall(($(prefix * "Cuhre"), libcuba), Cdouble,
              (Cint, # ndim
               Cint, # ncomp
               Ptr{Void}, # integrand
               Any, # userdata
               $CubaInt, # nvec
               Cdouble, # reltol
               Cdouble, # abstol
               Cint, # flags
               $CubaInt, # minevals
               $CubaInt, # maxevals
               Cint, # key
               Ptr{Cchar}, # statefile
               Ptr{Void}, # spin
               Ptr{Cint}, # nregions
               Ptr{$CubaInt}, # neval
               Ptr{Cint}, # fail
               Ptr{Cdouble}, # integral
               Ptr{Cdouble}, # error
               Ptr{Cdouble}),# prob
              # Input
              x.ndim, x.ncomp, integrand, x.func, x.nvec, x.reltol, x.abstol,
              x.flags, x.minevals, x.maxevals, x.key, x.statefile, x.spin,
              # Output
              nregions, neval, fail, integral, error, prob)
    end

    func = Symbol(prefix, "cuhre")
    @eval function $func{T}(integrand::T, ndim::Integer=1, ncomp::Integer=1;
                            nvec::Integer=NVEC, reltol::Real=RELTOL, abstol::Real=ABSTOL,
                            flags::Integer=FLAGS, minevals::Real=MINEVALS,
                            maxevals::Real=MAXEVALS, key::Integer=KEY,
                            statefile::AbstractString=STATEFILE, spin::Ptr{Void}=SPIN)
        # Cuhre requires "ndim" to be at least 2, even for an integral over a one
        # dimensional domain.  Instead, we don't prevent users from setting wrong
        # "ndim" values like 0 or negative ones.
        ndim == 1 && (ndim = 2)
        return dointegrate(Cuhre(integrand, ndim, ncomp, $CubaInt(nvec), Cdouble(reltol),
                                 Cdouble(abstol), flags, trunc($CubaInt, minevals),
                                 trunc($CubaInt, maxevals), key, String(statefile), spin))
    end
end

"""
    cuhre(integrand, ndim=1, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions
    llcuhre(integrand, ndim=1, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions

Calculate integral of `integrand` over the unit hypercube in `ndim` dimensions
using Cuhre algorithm.  `integrand` is a vectorial function with `ncomp`
components.  `ndim` and `ncomp` default to 1.

Accepted keywords:

* `nvec`
* `reltol`
* `abstol`
* `flags`
* `minevals`
* `maxevals`
* `key`
* `statefile`
* `spin`
"""
cuhre, llcuhre
