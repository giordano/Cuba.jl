### divonne.jl --- Integrate with Divonne method

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

struct Divonne{T, D} <: Integrand{T}
    func::T
    userdata::D
    ndim::Int
    ncomp::Int
    nvec::Int64
    rtol::Cdouble
    atol::Cdouble
    flags::Int
    seed::Int
    minevals::Int64
    maxevals::Int64
    key1::Int
    key2::Int
    key3::Int
    maxpass::Int
    border::Cdouble
    maxchisq::Cdouble
    mindeviation::Cdouble
    ngiven::Int64
    ldxgiven::Int
    xgiven::Array{Cdouble, 2}
    nextra::Int64
    peakfinder::Ptr{Cvoid}
    statefile::String
    spin::Ptr{Cvoid}
end

@inline function dointegrate!(x::Divonne{T, D}, integrand, integral,
                              error, prob, neval, fail, nregions) where {T, D}

    userdata = ismissing(x.userdata) ? x.func : (x.func, x.userdata)

    ccall((:llDivonne, libcuba), Cdouble,
          (Cint, # ndim
           Cint, # ncomp
           Ptr{Cvoid}, # integrand
           Any, # userdata
           Int64, # nvec
           Cdouble, # rtol
           Cdouble, # atol
           Cint, # flags
           Cint, # seed
           Int64, # minevals
           Int64, # maxevals
           Cint, # key1
           Cint, # key2
           Cint, # key3
           Cint, # maxpass
           Cdouble, # border
           Cdouble, # maxchisq
           Cdouble, # mindeviation
           Int64, # ngiven
           Cint, # ldxgiven
           Ptr{Cdouble}, # xgiven
           Int64, # nextra
           Ptr{Cvoid}, # peakfinder
           Ptr{Cchar}, # statefile
           Ptr{Cvoid}, # spin
           Ptr{Cint}, # nregions
           Ptr{Int64}, # neval
           Ptr{Cint}, # fail
           Ptr{Cdouble}, # integral
           Ptr{Cdouble}, # error
           Ptr{Cdouble}),# prob
          # Input
          x.ndim, x.ncomp, integrand, userdata, x.nvec, x.rtol,
          x.atol, x.flags, x.seed, x.minevals, x.maxevals, x.key1, x.key2,
          x.key3, x.maxpass, x.border, x.maxchisq, x.mindeviation, x.ngiven,
          x.ldxgiven, x.xgiven, x.nextra, x.peakfinder, x.statefile, x.spin,
          # Output
          nregions, neval, fail, integral, error, prob)
end

"""
    divonne(integrand, ndim=2, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions

Calculate integral of `integrand` over the unit hypercube in `ndim` dimensions
using Divonne algorithm.  `integrand` is a vectorial function with `ncomp`
components. `ncomp` defaults to 1, `ndim` defaults to 2 and must be ≥ 2.

Accepted keywords:

* `userdata`
* `nvec`
* `rtol`
* `atol`
* `flags`
* `seed`
* `minevals`
* `maxevals`
* `key1`
* `key2`
* `key3`
* `maxpass`
* `border`
* `maxchisq`
* `mindeviation`
* `ngiven`
* `ldxgiven`
* `xgiven`
* `nextra`
* `peakfinder`
* `statefile`
* `spin`
"""
function divonne(integrand::T, ndim::Integer=2, ncomp::Integer=1;
                 nvec::Integer=NVEC, rtol::Real=RTOL,
                 atol::Real=ATOL, flags::Integer=FLAGS,
                 seed::Integer=SEED, minevals::Real=MINEVALS,
                 maxevals::Real=MAXEVALS, key1::Integer=KEY1,
                 key2::Integer=KEY2, key3::Integer=KEY3,
                 maxpass::Integer=MAXPASS, border::Real=BORDER,
                 maxchisq::Real=MAXCHISQ,
                 mindeviation::Real=MINDEVIATION,
                 ngiven::Integer=NGIVEN, ldxgiven::Integer=LDXGIVEN,
                 xgiven::Array{Cdouble,2}=zeros(Cdouble, ldxgiven,
                                                ngiven),
                 nextra::Integer=NEXTRA,
                 peakfinder::Ptr{Cvoid}=PEAKFINDER,
                 statefile::AbstractString=STATEFILE,
                 spin::Ptr{Cvoid}=SPIN, reltol=missing, abstol=missing, userdata=missing) where {T}
    atol_,rtol_ = tols(atol,rtol,abstol,reltol)
    # Divonne requires "ndim" to be at least 2, even for an integral over a one
    # dimensional domain.  Instead, we don't prevent users from setting wrong
    # "ndim" values like 0 or negative ones.
    ndim >=2 || throw(ArgumentError("In Divonne ndim must be ≥ 2"))
    return dointegrate(Divonne(integrand, userdata, ndim, ncomp, Int64(nvec), Cdouble(rtol_),
                               Cdouble(atol_), flags, seed, trunc(Int64, minevals),
                               trunc(Int64, maxevals), key1, key2, key3, maxpass,
                               Cdouble(border), Cdouble(maxchisq), Cdouble(mindeviation),
                               Int64(ngiven), ldxgiven, xgiven, Int64(nextra),
                               peakfinder, String(statefile), spin))
end
