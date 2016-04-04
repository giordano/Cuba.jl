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

export Vegas

const libcuba = joinpath(Pkg.dir("Cuba"), "deps", "libcuba")

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
                   ndim, ncomp, integrand_ptr, userdata, nvec, epsrel, epsabs,
                   verbose, seed, mineval, maxeval, nstart, nincrease,
                   nbatch, gridno, statefile, spin, neval, fail,
                   integral, error, prob
                   )

    return integral, error, prob, neval[], fail[]
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
      ) = Vegas(integrand, ndim, ncomp, userdata, nvec, epsrel, epsabs,
                verbose, seed, trunc(Integer, mineval), trunc(Integer, maxeval),
                nstart, nincrease, nbatch, gridno, statefile, spin)

### Example of integrand function, just for test
function func(ndim::Cint,
              xx::Ptr{Cdouble},
              ncomp::Cint,
              ff::Ptr{Cdouble},
              userdata::Ptr{Void}=USERDATA_DEF
              )
    x = pointer_to_array(xx, (ndim,))
    f = pointer_to_array(ff, (ncomp,))
    for i = 1:ncomp
        f[i] = sin(x[1])*cos(x[2])*exp(x[3])
    end
    xx = pointer_from_objref(x)
    ff = pointer_from_objref(f)
    return Cint(0)::Cint
end

end # module
