### suave.jl --- Integrate with Suave method

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

struct Suave{T} <: Integrand{T}
    func::T
    ndim::Int
    ncomp::Int
    nvec::Int64
    reltol::Cdouble
    abstol::Cdouble
    flags::Int
    seed::Int
    minevals::Int64
    maxevals::Int64
    nnew::Int64
    nmin::Int64
    flatness::Cdouble
    statefile::String
    spin::Ptr{Cvoid}
end

@inline function dointegrate!(x::Suave{T}, integrand, integral,
                              error, prob, neval, fail, nregions) where {T}
    ccall((:llSuave, libcuba), Cdouble,
          (Cint, # ndim
           Cint, # ncomp
           Ptr{Cvoid}, # integrand
           Any, # userdata
           Int64, # nvec
           Cdouble, # reltol
           Cdouble, # abstol
           Cint, # flags
           Cint, # seed
           Int64, # minevals
           Int64, # maxevals
           Int64, # nnew
           Int64, # nmin
           Cdouble, # flatness
           Ptr{Cchar}, # statefile
           Ptr{Cvoid}, # spin
           Ptr{Cint}, # nregions
           Ptr{Int64}, # neval
           Ptr{Cint}, # fail
           Ptr{Cdouble}, # integral
           Ptr{Cdouble}, # error
           Ptr{Cdouble}),# prob
          # Input
          x.ndim, x.ncomp, integrand, x.func, x.nvec,
          x.reltol, x.abstol, x.flags, x.seed, x.minevals, x.maxevals,
          x.nnew, x.nmin, x.flatness, x.statefile, x.spin,
          # Output
          nregions, neval, fail, integral, error, prob)
end

"""
    suave(integrand, ndim=1, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions

Calculate integral of `integrand` over the unit hypercube in `ndim` dimensions
using Suave algorithm.  `integrand` is a vectorial function with `ncomp`
components. `ndim` and `ncomp` default to 1.

Accepted keywords:

* `nvec`
* `reltol`
* `abstol`
* `flags`
* `seed`
* `minevals`
* `maxevals`
* `nnew`
* `nmin`
* `flatness`
* `statefile`
* `spin`
"""
function suave(integrand::T, ndim::Integer=1, ncomp::Integer=1;
               nvec::Integer=NVEC, reltol::Real=RELTOL, abstol::Real=ABSTOL,
               flags::Integer=FLAGS, seed::Integer=SEED,
               minevals::Real=MINEVALS, maxevals::Real=MAXEVALS,
               nnew::Integer=NNEW, nmin::Integer=NMIN, flatness::Real=FLATNESS,
               statefile::AbstractString=STATEFILE, spin::Ptr{Cvoid}=SPIN) where {T}
    return dointegrate(Suave(integrand, ndim, ncomp, Int64(nvec), Cdouble(reltol),
                             Cdouble(abstol), flags, seed, trunc(Int64, minevals),
                             trunc(Int64, maxevals), Int64(nnew), Int64(nmin),
                             Cdouble(flatness), String(statefile), spin))
end
