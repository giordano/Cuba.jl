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
NDIM_DEF      = 3
NCOMP_DEF     = 1
USERDATA_DEF  = C_NULL
NVEC_DEF      = 1
EPSREL_DEF    = 1e-4
EPSABS_DEF    = 1e-12
VERBOSE_DEF   = 0
SEED_DEF      = 0
MINEVAL_DEF   = 0
MAXEVAL_DEF   = 1000000
NSTART_DEF    = 1000
NINCREASE_DEF = 500
NBATCH_DEF    = 1000
GRIDNO_DEF    = 0
STATEFILE_DEF = ""
SPIN_DEF      = C_NULL

NNEW_DEF     = 1000
NMIN_DEF     = 2
FLATNESS_DEF = 25.

KEY1_DEF         = 47
KEY2_DEF         = 1
KEY3_DEF         = 1
MAXPASS_DEF      = 5
BORDER_DEF       = 0.
MAXCHISQ_DEF     = 10.
MINDEVIATION_DEF = .25
NGIVEN_DEF       = 0
LDXGIVEN_DEF     = 0
XGIVEN_DEF       = zeros(typeof(1.0), LDXGIVEN_DEF, NGIVEN_DEF)
NEXTRA_DEF       = 0
PEAKFINDER_DEF   = C_NULL

KEY_DEF = 0

### Vegas
function Vegas(integrand::Function,
               ndim::Integer,
               ncomp::Integer,
               userdata::Ptr{Void},
               nvec::Integer,
               epsrel::Real,
               epsabs::Real,
               verbose::Integer,
               seed::Integer,
               mineval::Integer,
               maxeval::Integer,
               nstart::Integer,
               nincrease::Integer,
               nbatch::Integer,
               gridno::Integer,
               statefile::AbstractString,
               spin::Ptr{Void}
               )

    neval    = Ref{Cint}(0)
    fail     = Ref{Cint}(0)
    integral = zeros(typeof(1.0), ncomp)
    error    = zeros(typeof(1.0), ncomp)
    prob     = zeros(typeof(1.0), ncomp)

    const integrand_ptr = cfunction(integrand, Cint,
                                    (Ref{Cint}, # ndim
                                     Ptr{Cdouble}, # x
                                     Ref{Cint}, # ncomp
                                     Ptr{Cdouble}, # f
                                     Ptr{Void} # userdata
                                     ))

    result = ccall((:Vegas, libcuba), Cdouble,
                   (Cint, # ndim
                    Cint, # ncomp
                    Ptr{Void}, # integrand
                    Ptr{Void}, # userdata
                    Cint, # nvec
                    Cdouble, # epsrel
                    Cdouble, # epsabs
                    Cint, # verbose
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
                    Ptr{Cdouble}  # prob
                    ),
                   # Input
                   ndim, ncomp, integrand_ptr, userdata, nvec, epsrel, epsabs,
                   verbose, seed, mineval, maxeval, nstart, nincrease,
                   nbatch, gridno, statefile, spin,
                   # Output
                   neval, fail,
                   integral, error, prob
                   )

    return integral, error, prob, neval[], fail[], 0
end

Vegas(integrand::Function;
      ndim::Integer=NDIM_DEF,
      ncomp::Integer=NCOMP_DEF,
      userdata::Ptr{Void}=USERDATA_DEF,
      nvec::Integer=NVEC_DEF,
      epsrel::Real=EPSREL_DEF,
      epsabs::Real=EPSABS_DEF,
      verbose::Integer=VERBOSE_DEF,
      seed::Integer=SEED_DEF,
      mineval::Real=MINEVAL_DEF,
      maxeval::Real=MAXEVAL_DEF,
      nstart::Integer=NSTART_DEF,
      nincrease::Integer=NINCREASE_DEF,
      nbatch::Integer=NBATCH_DEF,
      gridno::Integer=GRIDNO_DEF,
      statefile::AbstractString=STATEFILE_DEF,
      spin::Ptr{Void}=SPIN_DEF
      ) = Vegas(integrand, ndim, ncomp, userdata, nvec, float(epsrel),
                float(epsabs),
                verbose, seed, trunc(Integer, mineval), trunc(Integer, maxeval),
                nstart, nincrease, nbatch, gridno, statefile, spin)

### Suave
function Suave(integrand::Function,
               ndim::Integer,
               ncomp::Integer,
               userdata::Ptr{Void},
               nvec::Integer,
               epsrel::Real,
               epsabs::Real,
               verbose::Integer,
               seed::Integer,
               mineval::Integer,
               maxeval::Integer,
               nnew::Integer,
               nmin::Integer,
               flatness::AbstractFloat,
               statefile::AbstractString,
               spin::Ptr{Void}
               )

    nregions = Ref{Cdouble}(0.0)
    neval    = Ref{Cint}(0)
    fail     = Ref{Cint}(0)
    integral = zeros(typeof(1.0), ncomp)
    error    = zeros(typeof(1.0), ncomp)
    prob     = zeros(typeof(1.0), ncomp)

    const integrand_ptr = cfunction(integrand, Cint,
                                    (Ref{Cint}, # ndim
                                     Ptr{Cdouble}, # x
                                     Ref{Cint}, # ncomp
                                     Ptr{Cdouble}, # f
                                     Ptr{Void} # userdata
                                     ))

    result = ccall((:Suave, libcuba), Cdouble,
                   (Cint, # ndim
                    Cint, # ncomp
                    Ptr{Void}, # integrand
                    Ptr{Void}, # userdata
                    Cint, # nvec
                    Cdouble, # epsrel
                    Cdouble, # epsabs
                    Cint, # verbose
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
                    Ptr{Cdouble}  # prob
                    ),
                   # Input
                   ndim, ncomp, integrand_ptr, userdata, nvec, epsrel, epsabs,
                   verbose, seed, mineval, maxeval, nnew, nmin,
                   flatness, statefile, spin,
                   # Output
                   nregions, neval, fail,
                   integral, error, prob
                   )

    return integral, error, prob, neval[], fail[], nregions[]
end

Suave(integrand::Function;
      ndim::Integer=NDIM_DEF,
      ncomp::Integer=NCOMP_DEF,
      userdata::Ptr{Void}=USERDATA_DEF,
      nvec::Integer=NVEC_DEF,
      epsrel::Real=EPSREL_DEF,
      epsabs::Real=EPSABS_DEF,
      verbose::Integer=VERBOSE_DEF,
      seed::Integer=SEED_DEF,
      mineval::Real=MINEVAL_DEF,
      maxeval::Real=MAXEVAL_DEF,
      nnew::Integer=NNEW_DEF,
      nmin::Integer=NMIN_DEF,
      flatness::Real=FLATNESS_DEF,
      statefile::AbstractString=STATEFILE_DEF,
      spin::Ptr{Void}=SPIN_DEF
      ) = Suave(integrand, ndim, ncomp, userdata, nvec, float(epsrel), float(epsabs),
                verbose, seed, trunc(Integer, mineval), trunc(Integer, maxeval),
                nnew, nmin, float(flatness), statefile, spin)

### Divonne
function Divonne{F<:AbstractFloat}(integrand::Function,
                                   ndim::Integer,
                                   ncomp::Integer,
                                   userdata::Ptr{Void},
                                   nvec::Integer,
                                   epsrel::Real,
                                   epsabs::Real,
                                   verbose::Integer,
                                   seed::Integer,
                                   mineval::Integer,
                                   maxeval::Integer,
                                   key1::Integer,
                                   key2::Integer,
                                   key3::Integer,
                                   maxpass::Integer,
                                   border::F,
                                   maxchisq::F,
                                   mindeviation::F,
                                   ngiven::Integer,
                                   ldxgiven::Integer,
                                   xgiven::AbstractArray{F},
                                   nextra::Integer,
                                   peakfinder::Ptr{Void},
                                   statefile::AbstractString,
                                   spin::Ptr{Void}
                                   )

    nregions = Ref{Cdouble}(0.0)
    neval    = Ref{Cint}(0)
    fail     = Ref{Cint}(0)
    integral = zeros(typeof(1.0), ncomp)
    error    = zeros(typeof(1.0), ncomp)
    prob     = zeros(typeof(1.0), ncomp)

    const integrand_ptr = cfunction(integrand, Cint,
                                    (Ref{Cint}, # ndim
                                     Ptr{Cdouble}, # x
                                     Ref{Cint}, # ncomp
                                     Ptr{Cdouble}, # f
                                     Ptr{Void} # userdata
                                     ))

    result = ccall((:Divonne, libcuba), Cdouble,
                   (Cint, # ndim
                    Cint, # ncomp
                    Ptr{Void}, # integrand
                    Ptr{Void}, # userdata
                    Cint, # nvec
                    Cdouble, # epsrel
                    Cdouble, # epsabs
                    Cint, # verbose
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
                    Ptr{Cdouble}  # prob
                    ),
                   # Input
                   ndim, ncomp, integrand_ptr, userdata, nvec, epsrel, epsabs,
                   verbose, seed, mineval, maxeval, key1, key2, key3, maxpass,
                   border, maxchisq, mindeviation, ngiven, ldxgiven, xgiven,
                   nextra, peakfinder, statefile, spin,
                   # Output
                   nregions, neval, fail,
                   integral, error, prob
                   )

    return integral, error, prob, neval[], fail[], nregions[]
end

Divonne{R<:Real}(integrand::Function;
                 ndim::Integer=NDIM_DEF,
                 ncomp::Integer=NCOMP_DEF,
                 userdata::Ptr{Void}=USERDATA_DEF,
                 nvec::Integer=NVEC_DEF,
                 epsrel::Real=EPSREL_DEF,
                 epsabs::Real=EPSABS_DEF,
                 verbose::Integer=VERBOSE_DEF,
                 seed::Integer=SEED_DEF,
                 mineval::Real=MINEVAL_DEF,
                 maxeval::Real=MAXEVAL_DEF,
                 key1::Integer=KEY1_DEF,
                 key2::Integer=KEY2_DEF,
                 key3::Integer=KEY3_DEF,
                 maxpass::Integer=MAXPASS_DEF,
                 border::Real=BORDER_DEF,
                 maxchisq::Real=MAXCHISQ_DEF,
                 mindeviation::Real=MINDEVIATION_DEF,
                 ngiven::Integer=NGIVEN_DEF,
                 ldxgiven::Integer=LDXGIVEN_DEF,
                 xgiven::AbstractArray{R}=XGIVEN_DEF,
                 nextra::Integer=NEXTRA_DEF,
                 peakfinder::Ptr{Void}=PEAKFINDER_DEF,
                 statefile::AbstractString=STATEFILE_DEF,
                 spin::Ptr{Void}=SPIN_DEF
                 ) = Divonne(integrand, ndim, ncomp, userdata, nvec,
                             float(epsrel), float(epsabs),
                             verbose, seed, trunc(Integer, mineval),
                             trunc(Integer, maxeval),
                             key1, key2, key3, maxpass,
                             float(border), float(maxchisq),
                             float(mindeviation), ngiven, ldxgiven,
                             xgiven,
                             nextra, peakfinder, statefile, spin)

### Cuhre
function Cuhre(integrand::Function,
               ndim::Integer,
               ncomp::Integer,
               userdata::Ptr{Void},
               nvec::Integer,
               epsrel::Real,
               epsabs::Real,
               verbose::Integer,
               mineval::Integer,
               maxeval::Integer,
               key::Integer,
               statefile::AbstractString,
               spin::Ptr{Void}
               )

    nregions = Ref{Cdouble}(0.0)
    neval    = Ref{Cint}(0)
    fail     = Ref{Cint}(0)
    integral = zeros(typeof(1.0), ncomp)
    error    = zeros(typeof(1.0), ncomp)
    prob     = zeros(typeof(1.0), ncomp)

    const integrand_ptr = cfunction(integrand, Cint,
                                    (Ref{Cint}, # ndim
                                     Ptr{Cdouble}, # x
                                     Ref{Cint}, # ncomp
                                     Ptr{Cdouble}, # f
                                     Ptr{Void} # userdata
                                     ))

    result = ccall((:Cuhre, libcuba), Cdouble,
                   (Cint, # ndim
                    Cint, # ncomp
                    Ptr{Void}, # integrand
                    Ptr{Void}, # userdata
                    Cint, # nvec
                    Cdouble, # epsrel
                    Cdouble, # epsabs
                    Cint, # verbose
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
                    Ptr{Cdouble}  # prob
                    ),
                   # Input
                   ndim, ncomp, integrand_ptr, userdata, nvec, epsrel, epsabs,
                   verbose, mineval, maxeval, key,
                   statefile, spin,
                   # Output
                   nregions, neval, fail,
                   integral, error, prob
                   )

    return integral, error, prob, neval[], fail[], nregions[]
end

Cuhre(integrand::Function;
      ndim::Integer=NDIM_DEF,
      ncomp::Integer=NCOMP_DEF,
      userdata::Ptr{Void}=USERDATA_DEF,
      nvec::Integer=NVEC_DEF,
      epsrel::Real=EPSREL_DEF,
      epsabs::Real=EPSABS_DEF,
      verbose::Integer=VERBOSE_DEF,
      mineval::Real=MINEVAL_DEF,
      maxeval::Real=MAXEVAL_DEF,
      key::Integer=KEY_DEF,
      statefile::AbstractString=STATEFILE_DEF,
      spin::Ptr{Void}=SPIN_DEF
      ) = Cuhre(integrand, ndim, ncomp, userdata, nvec, float(epsrel),
                float(epsabs),
                verbose, trunc(Integer, mineval), trunc(Integer, maxeval),
                key, statefile, spin)

end # module
