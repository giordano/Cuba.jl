### Cuba.jl --- Julia library for multidimensional numerical integration.

# Copyright (C) 2016  Mosè Giordano

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

__precompile__()

module Cuba

export vegas, suave, divonne, cuhre

# Note: don't use Pkg.dir("PkgName") here because the package may be installed
# elsewhere.
const libcuba = joinpath(dirname(@__FILE__), "..", "deps", "libcuba")

### Default values of parameters
# Common arguments.
const NVEC      = 1
const RELTOL    = 1e-4
const ABSTOL    = 1e-12
const FLAGS     = 0
const SEED      = 0
const MINEVALS  = 0
const MAXEVALS  = 1000000
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
const FLATNESS = 25.0

# Divonne-specific arguments.
const KEY1         = 47
const KEY2         = 1
const KEY3         = 1
const MAXPASS      = 5
const BORDER       = 0.0
const MAXCHISQ     = 10.0
const MINDEVIATION = 0.25
const NGIVEN       = 0
const LDXGIVEN     = 0
const XGIVEN       = 0
const NEXTRA       = 0
const PEAKFINDER   = C_NULL

# Cuhre-specific argument.
const KEY = 0

### Functions

# Return pointer for "integrand", to be passed as "integrand" argument to Cuba
# functions.
integrand_ptr(integrand::Function) = cfunction(integrand, Cint,
                                               (Ref{Cint}, # ndim
                                                Ptr{Cdouble}, # x
                                                Ref{Cint}, # ncomp
                                                Ptr{Cdouble}, # f
                                                Ptr{Void})) # userdata

# Note on implementation: instead of passing the function that performs
# calculations as "integrand" argument to integrator routines, we pass the
# pointer to this function and use "func_" to actually perform calculations.
# Without doing so and trying to "cfunction" a function not yet defined we would
# run into a
#
#   cfunction: no method exactly matched the required type signature (function not yet c-callable)
#
# error.  See http://julialang.org/blog/2013/05/callback for more information on
# this, in particular the section about "qsort_r" ("Passing closures via
# pass-through pointers").  Thanks to Steven G. Johnson for pointing to this.
function generic_integrand!(ndim::Cint, x_::Ptr{Cdouble}, ncomp::Cint,
                            f_::Ptr{Cdouble}, func_::Ptr{Void})
    # Get arrays from "x_" and "f_" pointers.
    x = unsafe_wrap(Array, x_, (ndim,))
    f = unsafe_wrap(Array, f_, (ncomp,))
    # Get the function from "func_" pointer.
    func! = unsafe_pointer_to_objref(func_)::Function
    func!(x, f)
    return Cint(0)
end

function __init__()
    # This pointer needs to be available at runtime.  See "Module initialization
    # and precompilation" section of Julia manual.
    global const c_generic_integrand! = integrand_ptr(generic_integrand!)
    Cuba.cores(0, 10000)
end

# One function to rule them all.
function dointegrate(algorithm::Symbol,
                     # First common arguments.
                     integrand::Ptr{Void}, ndim::Integer, ncomp::Integer,
                     userdata::Ptr{Void}, nvec::Integer, reltol::Real,
                     abstol::Real, flags::Integer, seed::Integer,
                     minevals::Integer, maxevals::Integer,
                     # Vegas-specific arguments.
                     nstart::Integer, nincrease::Integer,
                     nbatch::Integer, gridno::Integer,
                     # Suave-specific arguments.
                     nnew::Integer, nmin::Integer, flatness::Real,
                     # Divonne-specific arguments.
                     key1::Integer, key2::Integer, key3::Integer,
                     maxpass::Integer, border::Real, maxchisq::Real,
                     mindeviation::Real, ngiven::Integer, ldxgiven::Integer,
                     xgiven::Any, nextra::Integer, peakfinder::Ptr{Void},
                     # Cuhre-specific argument.
                     key::Integer,
                     # Final common arguments.
                     statefile::AbstractString, spin::Ptr{Void})
    nregions = Ref{Cint}(0)
    neval    = Ref{Cint}(0)
    fail     = Ref{Cint}(0)
    integral = zeros(Cdouble, ncomp)
    error    = zeros(Cdouble, ncomp)
    prob     = zeros(Cdouble, ncomp)
    if algorithm == :Cuhre
        ccall((:Cuhre, libcuba), Cdouble,
              (Cint, # ndim
               Cint, # ncomp
               Ptr{Void}, # integrand
               Ptr{Void}, # userdata
               Cint, # nvec
               Cdouble, # reltol
               Cdouble, # abstol
               Cint, # flags
               Cint, # minevals
               Cint, # maxevals
               Cint, # key
               Ptr{Cchar}, # statefile
               Ptr{Void}, # spin
               Ptr{Cint}, # nregions
               Ptr{Cint}, # neval
               Ptr{Cint}, # fail
               Ptr{Cdouble}, # integral
               Ptr{Cdouble}, # error
               Ptr{Cdouble}),# prob
              # Input
              ndim, ncomp, integrand, userdata, nvec, reltol,
              abstol, flags, minevals, maxevals, key, statefile, spin,
              # Output
              nregions, neval, fail, integral, error, prob)
    elseif algorithm == :Vegas
        ccall((:Vegas, libcuba), Cdouble,
              (Cint, # ndim
               Cint, # ncomp
               Ptr{Void}, # integrand
               Ptr{Void}, # userdata
               Cint, # nvec
               Cdouble, # reltol
               Cdouble, # abstol
               Cint, # flags
               Cint, # seed
               Cint, # minevals
               Cint, # maxevals
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
              ndim, ncomp, integrand, userdata, nvec,
              reltol, abstol, flags, seed, minevals, maxevals,
              nstart, nincrease, nbatch, gridno, statefile, spin,
              # Output
              neval, fail, integral, error, prob)
    elseif algorithm == :Suave
        ccall((:Suave, libcuba), Cdouble,
              (Cint, # ndim
               Cint, # ncomp
               Ptr{Void}, # integrand
               Ptr{Void}, # userdata
               Cint, # nvec
               Cdouble, # reltol
               Cdouble, # abstol
               Cint, # flags
               Cint, # seed
               Cint, # minevals
               Cint, # maxevals
               Cint, # nnew
               Cint, # nmin
               Cdouble, # flatness
               Ptr{Cchar}, # statefile
               Ptr{Void}, # spin
               Ptr{Cint}, # nregions
               Ptr{Cint}, # neval
               Ptr{Cint}, # fail
               Ptr{Cdouble}, # integral
               Ptr{Cdouble}, # error
               Ptr{Cdouble}),# prob
              # Input
              ndim, ncomp, integrand, userdata, nvec,
              reltol, abstol, flags, seed, minevals, maxevals,
              nnew, nmin, flatness, statefile, spin,
              # Output
              nregions, neval, fail, integral, error, prob)
    elseif algorithm == :Divonne
        ccall((:Divonne, libcuba), Cdouble,
              (Cint, # ndim
               Cint, # ncomp
               Ptr{Void}, # integrand
               Ptr{Void}, # userdata
               Cint, # nvec
               Cdouble, # reltol
               Cdouble, # abstol
               Cint, # flags
               Cint, # seed
               Cint, # minevals
               Cint, # maxevals
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
               Ptr{Cint}, # nregions
               Ptr{Cint}, # neval
               Ptr{Cint}, # fail
               Ptr{Cdouble}, # integral
               Ptr{Cdouble}, # error
               Ptr{Cdouble}),# prob
              # Input
              ndim, ncomp, integrand, userdata, nvec, reltol,
              abstol, flags, seed, minevals, maxevals, key1, key2, key3,
              maxpass, border, maxchisq, mindeviation, ngiven, ldxgiven,
              xgiven, nextra, peakfinder, statefile, spin,
              # Output
              nregions, neval, fail, integral, error, prob)
    end
    return integral, error, prob, neval[], fail[], nregions[]
end

"""
    vegas(integrand, ndim=1, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions

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
vegas(integrand::Function, ndim::Integer=1, ncomp::Integer=1;
      nvec::Integer=NVEC,
      reltol::Real=RELTOL, abstol::Real=ABSTOL, flags::Integer=FLAGS,
      seed::Integer=SEED, minevals::Real=MINEVALS, maxevals::Real=MAXEVALS,
      nstart::Integer=NSTART, nincrease::Integer=NINCREASE,
      nbatch::Integer=NBATCH, gridno::Integer=GRIDNO,
      statefile::AbstractString=STATEFILE, spin::Ptr{Void}=SPIN) =
          dointegrate(:Vegas, c_generic_integrand!::Ptr{Void}, ndim, ncomp,
                      pointer_from_objref(integrand), nvec, float(reltol),
                      float(abstol), flags, seed, trunc(Integer, minevals),
                      trunc(Integer, maxevals), nstart, nincrease, nbatch,
                      gridno, NNEW, NMIN, FLATNESS, KEY1, KEY2, KEY3, MAXPASS,
                      BORDER, MAXCHISQ, MINDEVIATION, NGIVEN, LDXGIVEN, XGIVEN,
                      NEXTRA, PEAKFINDER, KEY, statefile, spin)

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
suave(integrand::Function, ndim::Integer=1, ncomp::Integer=1;
      nvec::Integer=NVEC, reltol::Real=RELTOL,
      abstol::Real=ABSTOL, flags::Integer=FLAGS, seed::Integer=SEED,
      minevals::Real=MINEVALS, maxevals::Real=MAXEVALS, nnew::Integer=NNEW,
      nmin::Integer=NMIN, flatness::Real=FLATNESS,
      statefile::AbstractString=STATEFILE, spin::Ptr{Void}=SPIN) =
          dointegrate(:Suave, c_generic_integrand!::Ptr{Void}, ndim, ncomp,
                      pointer_from_objref(integrand), nvec, float(reltol),
                      float(abstol), flags, seed, trunc(Integer, minevals),
                      trunc(Integer, maxevals), NSTART, NINCREASE, NBATCH,
                      GRIDNO, nnew, nmin, flatness, KEY1, KEY2, KEY3, MAXPASS,
                      BORDER, MAXCHISQ, MINDEVIATION, NGIVEN, LDXGIVEN, XGIVEN,
                      NEXTRA, PEAKFINDER, KEY, statefile, spin)

"""
    divonne(integrand, ndim=1, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions

Calculate integral of `integrand` over the unit hypercube in `ndim` dimensions
using Divonne algorithm.  `integrand` is a vectorial function with `ncomp`
components. `ndim` and `ncomp` default to 1.

Accepted keywords:

* `nvec`
* `reltol`
* `abstol`
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
function divonne{R<:Real}(integrand::Function, ndim::Integer=1, ncomp::Integer=1;
                          nvec::Integer=NVEC, reltol::Real=RELTOL,
                          abstol::Real=ABSTOL, flags::Integer=FLAGS,
                          seed::Integer=SEED, minevals::Real=MINEVALS,
                          maxevals::Real=MAXEVALS, key1::Integer=KEY1,
                          key2::Integer=KEY2, key3::Integer=KEY3,
                          maxpass::Integer=MAXPASS, border::Real=BORDER,
                          maxchisq::Real=MAXCHISQ,
                          mindeviation::Real=MINDEVIATION,
                          ngiven::Integer=NGIVEN, ldxgiven::Integer=LDXGIVEN,
                          xgiven::AbstractArray{R}=zeros(Cdouble, ldxgiven,
                                                         ngiven),
                          nextra::Integer=NEXTRA,
                          peakfinder::Ptr{Void}=PEAKFINDER,
                          statefile::AbstractString=STATEFILE,
                          spin::Ptr{Void}=SPIN)
    # Divonne requires "ndim" to be at least 2, even for an integral over a one
    # dimensional domain.  Instead, we don't prevent users from setting wrong
    # "ndim" values like 0 or negative ones.
    ndim == 1 && (ndim = 2)
    return dointegrate(:Divonne, c_generic_integrand!::Ptr{Void}, ndim, ncomp,
                       pointer_from_objref(integrand), nvec, float(reltol),
                       float(abstol), flags, seed, trunc(Integer, minevals),
                       trunc(Integer, maxevals), NSTART, NINCREASE, NBATCH,
                       GRIDNO, NNEW, NMIN, FLATNESS, key1, key2, key3, maxpass,
                       float(border), float(maxchisq), float(mindeviation),
                       ngiven, ldxgiven, xgiven, nextra, peakfinder, KEY,
                       statefile, spin)
end

"""
    cuhre(integrand, ndim=1, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions

Calculate integral of `integrand` over the unit hypercube in `ndim` dimensions
using Cuhre algorithm.  `integrand` is a vectorial function with `ncomp`
components. `ndim` and `ncomp` default to 1.

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
function cuhre(integrand::Function, ndim::Integer=1, ncomp::Integer=1;
               nvec::Integer=NVEC, reltol::Real=RELTOL, abstol::Real=ABSTOL,
               flags::Integer=FLAGS, minevals::Real=MINEVALS,
               maxevals::Real=MAXEVALS, key::Integer=KEY,
               statefile::AbstractString=STATEFILE, spin::Ptr{Void}=SPIN)
    # Cuhre requires "ndim" to be at least 2, even for an integral over a one
    # dimensional domain.  Instead, we don't prevent users from setting wrong
    # "ndim" values like 0 or negative ones.
    ndim == 1 && (ndim = 2)
    return dointegrate(:Cuhre, c_generic_integrand!::Ptr{Void}, ndim, ncomp,
                       pointer_from_objref(integrand), nvec, float(reltol),
                       float(abstol), flags, SEED, trunc(Integer, minevals),
                       trunc(Integer, maxevals), NSTART, NINCREASE, NBATCH,
                       GRIDNO, NNEW, NMIN, FLATNESS, KEY1, KEY2, KEY3, MAXPASS,
                       BORDER, MAXCHISQ, MINDEVIATION, NGIVEN, LDXGIVEN, XGIVEN,
                       NEXTRA, PEAKFINDER, key, statefile, spin)
end

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
