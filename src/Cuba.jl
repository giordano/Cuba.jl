### Cuba.jl --- Julia library for multidimensional numerical integration.

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

module Cuba

export Vegas, Suave, Divonne, Cuhre

# Note: don't use Pkg.dir("PkgName") here because the package may be installed
# elsewhere.
const libcuba = joinpath(dirname(@__FILE__), "..", "deps", "libcuba")

### Default values of parameters
# Common arguments.
const USERDATA  = C_NULL
const NVEC      = 1
const EPSREL    = 1e-4
const EPSABS    = 1e-12
const FLAGS     = 0
const SEED      = 0
const MINEVAL   = 0
const MAXEVAL   = 1000000
const STATEFILE = ""
const SPIN      = C_NULL

# Vegas-specific arguments.
const NSTART    = 1000
const NINCREASE = 500
const NBATCH    = 1000
const GRIDNO    = 0

# Suave-specific arguments.
const NNEW     = 1000
const NMIN     = 2
const FLATNESS = 25.

# Divonne-specific arguments.
const KEY1         = 47
const KEY2         = 1
const KEY3         = 1
const MAXPASS      = 5
const BORDER       = 0.
const MAXCHISQ     = 10.
const MINDEVIATION = .25
const NGIVEN       = 0
const LDXGIVEN     = 0
const NEXTRA       = 0
const PEAKFINDER   = C_NULL

# Cuhre-specific argument.
const KEY = 0

### Functions

# Return pointer for "func", to be passed as "integrand" argument to Cuba
# functions.
function func_ptr(func::Function)
    return cfunction(func, Cint,
                     (Ref{Cint}, # ndim
                      Ptr{Cdouble}, # x
                      Ref{Cint}, # ncomp
                      Ptr{Cdouble}, # f
                      Ptr{Void})) # userdata
end

### Vegas
function Vegas(integrand::Function, ndim::Integer, ncomp::Integer,
               userdata::Ptr{Void}, nvec::Integer, epsrel::Real,
               epsabs::Real, flags::Integer, seed::Integer,
               mineval::Integer, maxeval::Integer, nstart::Integer,
               nincrease::Integer, nbatch::Integer, gridno::Integer,
               statefile::AbstractString, spin::Ptr{Void})
    Cuba.cores(0, 10000)
    neval    = Ref{Cint}(0)
    fail     = Ref{Cint}(0)
    integral = zeros(typeof(1.0), ncomp)
    error    = zeros(typeof(1.0), ncomp)
    prob     = zeros(typeof(1.0), ncomp)
    ccall((:Vegas, libcuba), Cdouble,
          (Cint, # ndim
           Cint, # ncomp
           Ptr{Void}, # integrand
           Ptr{Void}, # userdata
           Cint, # nvec
           Cdouble, # epsrel
           Cdouble, # epsabs
           Cint, # flags
           Cint, # seed
           Cint, # mineval
           Cint, # maxeval
           Cint, # nstart
           Cint, # nincrease
           Cint, # nbatch
           Cint, # gridno
           Ptr{Cchar}, # statefile
           Ptr{Void}, # spin
           Ptr{Cint}, # neval
           Ptr{Cint}, # fail
           Ptr{Cdouble}, # integral
           Ptr{Cdouble}, # error
           Ptr{Cdouble}),# prob
          # Input
          ndim, ncomp, func_ptr(integrand), userdata, nvec,
          epsrel, epsabs, flags, seed, mineval, maxeval,
          nstart, nincrease, nbatch, gridno, statefile, spin,
          # Output
          neval, fail, integral, error, prob)
    return integral, error, prob, neval[], fail[], 0
end

Vegas(integrand::Function, ndim::Integer, ncomp::Integer;
      userdata::Ptr{Void}=USERDATA, nvec::Integer=NVEC, epsrel::Real=EPSREL,
      epsabs::Real=EPSABS, flags::Integer=FLAGS, seed::Integer=SEED,
      mineval::Real=MINEVAL, maxeval::Real=MAXEVAL, nstart::Integer=NSTART,
      nincrease::Integer=NINCREASE, nbatch::Integer=NBATCH,
      gridno::Integer=GRIDNO, statefile::AbstractString=STATEFILE,
      spin::Ptr{Void}=SPIN) =
          Vegas(integrand, ndim, ncomp, userdata, nvec, float(epsrel),
                float(epsabs), flags, seed, trunc(Integer, mineval),
                trunc(Integer, maxeval), nstart, nincrease, nbatch,
                gridno, statefile, spin)

### Suave
function Suave(integrand::Function, ndim::Integer, ncomp::Integer,
               userdata::Ptr{Void}, nvec::Integer, epsrel::Real,
               epsabs::Real, flags::Integer, seed::Integer,
               mineval::Integer, maxeval::Integer, nnew::Integer,
               nmin::Integer, flatness::AbstractFloat,
               statefile::AbstractString, spin::Ptr{Void})
    Cuba.cores(0, 10000)
    nregions = Ref{Cdouble}(0.0)
    neval    = Ref{Cint}(0)
    fail     = Ref{Cint}(0)
    integral = zeros(typeof(1.0), ncomp)
    error    = zeros(typeof(1.0), ncomp)
    prob     = zeros(typeof(1.0), ncomp)
    ccall((:Suave, libcuba), Cdouble,
          (Cint, # ndim
           Cint, # ncomp
           Ptr{Void}, # integrand
           Ptr{Void}, # userdata
           Cint, # nvec
           Cdouble, # epsrel
           Cdouble, # epsabs
           Cint, # flags
           Cint, # seed
           Cint, # mineval
           Cint, # maxeval
           Cint, # nnew
           Cint, # nmin
           Cdouble, # flatness
           Ptr{Cchar}, # statefile
           Ptr{Void}, # spin
           Ptr{Cdouble}, # nregions
           Ptr{Cint}, # neval
           Ptr{Cint}, # fail
           Ptr{Cdouble}, # integral
           Ptr{Cdouble}, # error
           Ptr{Cdouble}),# prob
          # Input
          ndim, ncomp, func_ptr(integrand), userdata, nvec,
          epsrel, epsabs, flags, seed, mineval, maxeval,
          nnew, nmin, flatness, statefile, spin,
          # Output
          nregions, neval, fail, integral, error, prob)
    return integral, error, prob, neval[], fail[], nregions[]
end

Suave(integrand::Function, ndim::Integer, ncomp::Integer;
      userdata::Ptr{Void}=USERDATA, nvec::Integer=NVEC, epsrel::Real=EPSREL,
      epsabs::Real=EPSABS, flags::Integer=FLAGS, seed::Integer=SEED,
      mineval::Real=MINEVAL, maxeval::Real=MAXEVAL, nnew::Integer=NNEW,
      nmin::Integer=NMIN, flatness::Real=FLATNESS,
      statefile::AbstractString=STATEFILE, spin::Ptr{Void}=SPIN) =
          Suave(integrand, ndim, ncomp, userdata, nvec, float(epsrel),
                float(epsabs), flags, seed, trunc(Integer, mineval),
                trunc(Integer, maxeval), nnew, nmin, float(flatness),
                statefile, spin)

### Divonne
function Divonne{F<:AbstractFloat}(integrand::Function, ndim::Integer,
                                   ncomp::Integer, userdata::Ptr{Void},
                                   nvec::Integer, epsrel::Real, epsabs::Real,
                                   flags::Integer, seed::Integer,
                                   mineval::Integer, maxeval::Integer,
                                   key1::Integer, key2::Integer, key3::Integer,
                                   maxpass::Integer, border::F, maxchisq::F,
                                   mindeviation::F, ngiven::Integer,
                                   ldxgiven::Integer, xgiven::AbstractArray{F},
                                   nextra::Integer, peakfinder::Ptr{Void},
                                   statefile::AbstractString, spin::Ptr{Void})
    Cuba.cores(0, 10000)
    # Divonne requires "ndim" to be at least 2, even for an integral over a one
    # dimensional domain.  Instead, we don't prevent users from setting wrong
    # "ndim" values like 0 or negative ones.
    ndim == 1 && (ndim = 2)
    nregions = Ref{Cdouble}(0.0)
    neval    = Ref{Cint}(0)
    fail     = Ref{Cint}(0)
    integral = zeros(typeof(1.0), ncomp)
    error    = zeros(typeof(1.0), ncomp)
    prob     = zeros(typeof(1.0), ncomp)
    ccall((:Divonne, libcuba), Cdouble,
          (Cint, # ndim
           Cint, # ncomp
           Ptr{Void}, # integrand
           Ptr{Void}, # userdata
           Cint, # nvec
           Cdouble, # epsrel
           Cdouble, # epsabs
           Cint, # flags
           Cint, # seed
           Cint, # mineval
           Cint, # maxeval
           Cint, # key1
           Cint, # key2
           Cint, # key3
           Cint, # maxpass
           Cdouble, # border
           Cdouble, # maxchisq
           Cdouble, # mindeviation
           Cint, # ngiven
           Cint, # ldxgiven
           Ptr{Cdouble}, # xgiven
           Cint, # nextra
           Ptr{Void}, # peakfinder
           Ptr{Cchar}, # statefile
           Ptr{Void}, # spin
           Ptr{Cdouble}, # nregions
           Ptr{Cint}, # neval
           Ptr{Cint}, # fail
           Ptr{Cdouble}, # integral
           Ptr{Cdouble}, # error
           Ptr{Cdouble}),# prob
          # Input
          ndim, ncomp, func_ptr(integrand), userdata, nvec, epsrel,
          epsabs, flags, seed, mineval, maxeval, key1, key2, key3,
          maxpass, border, maxchisq, mindeviation, ngiven, ldxgiven,
          xgiven, nextra, peakfinder, statefile, spin,
          # Output
          nregions, neval, fail, integral, error, prob)
    return integral, error, prob, neval[], fail[], nregions[]
end

Divonne{R<:Real}(integrand::Function, ndim::Integer, ncomp::Integer;
                 userdata::Ptr{Void}=USERDATA, nvec::Integer=NVEC,
                 epsrel::Real=EPSREL, epsabs::Real=EPSABS,
                 flags::Integer=FLAGS, seed::Integer=SEED,
                 mineval::Real=MINEVAL, maxeval::Real=MAXEVAL,
                 key1::Integer=KEY1, key2::Integer=KEY2, key3::Integer=KEY3,
                 maxpass::Integer=MAXPASS, border::Real=BORDER,
                 maxchisq::Real=MAXCHISQ, mindeviation::Real=MINDEVIATION,
                 ngiven::Integer=NGIVEN, ldxgiven::Integer=LDXGIVEN,
                 xgiven::AbstractArray{R}=zeros(typeof(0.0), ldxgiven, ngiven),
                 nextra::Integer=NEXTRA, peakfinder::Ptr{Void}=PEAKFINDER,
                 statefile::AbstractString=STATEFILE, spin::Ptr{Void}=SPIN) =
                     Divonne(integrand, ndim, ncomp, userdata, nvec,
                             float(epsrel), float(epsabs), flags, seed,
                             trunc(Integer, mineval), trunc(Integer, maxeval),
                             key1, key2, key3, maxpass, float(border),
                             float(maxchisq), float(mindeviation), ngiven,
                             ldxgiven, xgiven, nextra, peakfinder, statefile,
                             spin)

### Cuhre
function Cuhre(integrand::Function, ndim::Integer, ncomp::Integer,
               userdata::Ptr{Void}, nvec::Integer, epsrel::Real, epsabs::Real,
               flags::Integer, mineval::Integer, maxeval::Integer,
               key::Integer, statefile::AbstractString, spin::Ptr{Void})
    Cuba.cores(0, 10000)
    # Cuhre requires "ndim" to be at least 2, even for an integral over a one
    # dimensional domain.  Instead, we don't prevent users from setting wrong
    # "ndim" values like 0 or negative ones.
    ndim == 1 && (ndim = 2)
    nregions = Ref{Cdouble}(0.0)
    neval    = Ref{Cint}(0)
    fail     = Ref{Cint}(0)
    integral = zeros(typeof(1.0), ncomp)
    error    = zeros(typeof(1.0), ncomp)
    prob     = zeros(typeof(1.0), ncomp)
    ccall((:Cuhre, libcuba), Cdouble,
          (Cint, # ndim
           Cint, # ncomp
           Ptr{Void}, # integrand
           Ptr{Void}, # userdata
           Cint, # nvec
           Cdouble, # epsrel
           Cdouble, # epsabs
           Cint, # flags
           Cint, # mineval
           Cint, # maxeval
           Cint, # key
           Ptr{Cchar}, # statefile
           Ptr{Void}, # spin
           Ptr{Cdouble}, # nregions
           Ptr{Cint}, # neval
           Ptr{Cint}, # fail
           Ptr{Cdouble}, # integral
           Ptr{Cdouble}, # error
           Ptr{Cdouble}),# prob
          # Input
          ndim, ncomp, func_ptr(integrand), userdata, nvec, epsrel,
          epsabs, flags, mineval, maxeval, key, statefile, spin,
          # Output
          nregions, neval, fail, integral, error, prob)
    return integral, error, prob, neval[], fail[], nregions[]
end

Cuhre(integrand::Function, ndim::Integer, ncomp::Integer;
      userdata::Ptr{Void}=USERDATA, nvec::Integer=NVEC, epsrel::Real=EPSREL,
      epsabs::Real=EPSABS, flags::Integer=FLAGS, mineval::Real=MINEVAL,
      maxeval::Real=MAXEVAL, key::Integer=KEY,
      statefile::AbstractString=STATEFILE, spin::Ptr{Void}=SPIN) =
          Cuhre(integrand, ndim, ncomp, userdata, nvec, float(epsrel),
                float(epsabs), flags, trunc(Integer, mineval),
                trunc(Integer, maxeval), key, statefile, spin)

### Other functions, not exported
function cores(n::Integer, p::Integer)
    ccall((:cubacores, libcuba), Ptr{Void}, (Cint, Cint), n, p)
    return 0
end

function accel(n::Integer, p::Integer)
    ccall((:cubaaccel, libcuba), Ptr{Void}, (Cint, Cint), n, p)
    return 0
end

function init(f::Ptr{Void}, arg::Ptr{Void})
    ccall((:cubainit, libcuba), Ptr{Void}, (Ptr{Void}, Ptr{Void}), f, arg)
    return 0
end

function exit(f::Ptr{Void}, arg::Ptr{Void})
    ccall((:cubaexit, libcuba), Ptr{Void}, (Ptr{Void}, Ptr{Void}), f, arg)
    return 0
end

end # module
