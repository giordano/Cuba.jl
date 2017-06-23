### Cuba.jl --- Julia library for multidimensional numerical integration.

# Copyright (C) 2016, 2017  Mosè Giordano

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
                            f_::Ptr{Cdouble}, func!)
    # Get arrays from "x_" and "f_" pointers.
    x = unsafe_wrap(Array, x_, (ndim,))
    f = unsafe_wrap(Array, f_, (ncomp,))
    func!(x, f)
    return Cint(0)
end

# Return pointer for "integrand", to be passed as "integrand" argument to Cuba functions.
integrand_ptr{T}(integrand::T) = cfunction(generic_integrand!, Cint,
                                           (Ref{Cint}, # ndim
                                            Ptr{Cdouble}, # x
                                            Ref{Cint}, # ncomp
                                            Ptr{Cdouble}, # f
                                            Ref{T})) # userdata

abstract Integrand{T}

function __init__()
    Cuba.cores(0, 10000)
end

include("cuhre.jl")
include("divonne.jl")
include("suave.jl")
include("vegas.jl")

immutable Integral
    integral::Vector{Float64}
    error::Vector{Float64}
    probability::Vector{Float64}
    neval::Int64
    fail::Int32
    nregions::Int32
end

Base.getindex(x::Integral, n::Integer) = getfield(x, n)
Base.start(x::Integral)   = 1
Base.next(x::Integral, i) = (x[i], i + 1)
Base.done(x::Integral, i) = (i > 6)

function Base.show(io::IO, x::Integral)
    ncomp = length(x.integral)
    println(io, ncomp == 1 ? "Component:" : "Components:")
    for i in 1:ncomp
        println(io, " ", lpad("$i", ceil(Int, log10(ncomp+1))), ": ", x.integral[i],
                " ± ", x.error[i], " (prob.: ", x.probability[i],")")
    end
    println(io, "Integrand evaluations: ", x.neval)
    println(io, "Fail:                  ", x.fail)
    print(io, "Number of subregions:  ", x.nregions)
end

@inline function dointegrate{T}(x::Integrand{T})
    integrand = integrand_ptr(x.func)
    integral  = Vector{Cdouble}(x.ncomp)
    error     = Vector{Cdouble}(x.ncomp)
    prob      = Vector{Cdouble}(x.ncomp)
    neval     = Ref{Int64}(0)
    fail      = Ref{Cint}(0)
    nregions  = Ref{Cint}(0)
    dointegrate!(x, integrand, integral, error, prob, neval, fail, nregions)
    return Integral(integral, error, prob, neval[], fail[], nregions[])
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

for f in (:vegas, :suave, :divonne, :cuhre)
    llf = Symbol("ll", f)
    @eval @deprecate $llf(args...; kwargs...) $f(args...; kwargs...)
end

end # module
