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

struct Cuhre{T} <: Integrand{T}
    func::T
    ndim::Int
    ncomp::Int
    nvec::Int64
    reltol::Cdouble
    abstol::Cdouble
    flags::Int
    minevals::Int64
    maxevals::Int64
    key::Int
    statefile::String
    spin::Ptr{Cvoid}
end

@inline function dointegrate!(x::Cuhre{T}, integrand, integral,
                              error, prob, neval, fail, nregions) where {T}
    ccall((:llCuhre, libcuba), Cdouble,
          (Cint, # ndim
           Cint, # ncomp
           Ptr{Cvoid}, # integrand
           Any, # userdata
           Int64, # nvec
           Cdouble, # reltol
           Cdouble, # abstol
           Cint, # flags
           Int64, # minevals
           Int64, # maxevals
           Cint, # key
           Ptr{Cchar}, # statefile
           Ptr{Cvoid}, # spin
           Ptr{Cint}, # nregions
           Ptr{Int64}, # neval
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

"""
    cuhre(integrand, ndim=2, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions

Calculate integral of `integrand` over the unit hypercube in `ndim` dimensions
using Cuhre algorithm.  `integrand` is a vectorial function with `ncomp`
components.  `ncomp` defaults to 1, `ndim` defaults to 2 and must be ≥ 2.

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
function cuhre(integrand::T, ndim::Integer=2, ncomp::Integer=1;
               nvec::Integer=NVEC, reltol::Real=RELTOL, abstol::Real=ABSTOL,
               flags::Integer=FLAGS, minevals::Real=MINEVALS,
               maxevals::Real=MAXEVALS, key::Integer=KEY,
               statefile::AbstractString=STATEFILE, spin::Ptr{Cvoid}=SPIN) where {T}
    # Cuhre requires "ndim" to be at least 2, even for an integral over a one
    # dimensional domain.  Instead, we don't prevent users from setting wrong
    # "ndim" values like 0 or negative ones.
    ndim >=2 || throw(ArgumentError("In Cuhre ndim must be ≥ 2"))
    return dointegrate(Cuhre(integrand, ndim, ncomp, Int64(nvec), Cdouble(reltol),
                             Cdouble(abstol), flags, trunc(Int64, minevals),
                             trunc(Int64, maxevals), key, String(statefile), spin))
end
