### Cuba.jl --- Julia library for multidimensional numerical integration.

__precompile__()

module Cuba

using Cuba_jll

export vegas, suave, divonne, cuhre

### Default values of parameters
# Common arguments.
const NVEC      = 1
const RTOL      = 1e-4
const ATOL      = 1e-12
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
#
# For better performance, when nvec == 1 we store a simple Vector, instead of a Matrix with
# second dimension equal to 1.
function generic_integrand!(ndim::Cint, x_::Ptr{Cdouble}, ncomp::Cint,
                            f_::Ptr{Cdouble}, func!)
    # Get arrays from "x_" and "f_" pointers.
    x = unsafe_wrap(Array, x_, (ndim,))
    f = unsafe_wrap(Array, f_, (ncomp,))
    func!(x, f)
    return Cint(0)
end
function generic_integrand!(ndim::Cint, x_::Ptr{Cdouble}, ncomp::Cint,
                            f_::Ptr{Cdouble}, func!, nvec::Cint)
    # Get arrays from "x_" and "f_" pointers.
    x = unsafe_wrap(Array, x_, (ndim, nvec))
    f = unsafe_wrap(Array, f_, (ncomp, nvec))
    func!(x, f)
    return Cint(0)
end

# Return pointer for "integrand", to be passed as "integrand" argument to Cuba functions.
integrand_ptr(integrand::T)  where {T} = @cfunction(generic_integrand!, Cint,
                                                    (Ref{Cint}, # ndim
                                                     Ptr{Cdouble}, # x
                                                     Ref{Cint}, # ncomp
                                                     Ptr{Cdouble}, # f
                                                     Ref{T})) # userdata
integrand_ptr_nvec(integrand::T) where {T} = @cfunction(generic_integrand!, Cint,
                                                        (Ref{Cint}, # ndim
                                                         Ptr{Cdouble}, # x
                                                         Ref{Cint}, # ncomp
                                                         Ptr{Cdouble}, # f
                                                         Ref{T}, # userdata
                                                         Ref{Cint})) # nvec

abstract type Integrand{T} end

function __init__()
    Cuba.cores(0, 10000)
end

# handle keyword deprecation
function tols(atol,rtol,abstol,reltol)
    if !ismissing(abstol) || !ismissing(reltol)
        Base.depwarn("abstol and reltol keywords are now atol and rtol, respectively", :quadgk)
    end
    return coalesce(abstol,atol), coalesce(reltol,rtol)
end

include("cuhre.jl")
include("divonne.jl")
include("suave.jl")
include("vegas.jl")

struct Integral{T}
    integral::Vector{Float64}
    error::Vector{Float64}
    probability::Vector{Float64}
    neval::Int64
    fail::Int32
    nregions::Int32
end

Base.getindex(x::Integral, n::Integer) = getfield(x, n)
function Base.iterate(x::Integral, i=1)
    i > 6 && return nothing
    return x[i], i + 1
end

_print_fail_extra(io::IO, x::Integral) = nothing
_print_fail_extra(io::IO, x::Integral{<:Divonne}) =
    print(io, "\nHint: Try increasing `maxevals` to ", x.neval+x.fail)
function print_fail(io::IO, x::Integral)
    print(io, "Note: ")
    fail = x.fail
    if fail < 0
        print(io, "Dimension out of range")
    elseif fail == 0
        print(io, "The desired accuracy was reached")
    elseif fail > 0
        print(io, "The accuracy was not met within the maximum number of evaluations")
        _print_fail_extra(io, x)
    end
end

function Base.show(io::IO, x::Integral)
    ncomp = length(x.integral)
    println(io, ncomp == 1 ? "Component:" : "Components:")
    for i in 1:ncomp
        println(io, " ", lpad("$i", ceil(Int, log10(ncomp+1))), ": ", x.integral[i],
                " ± ", x.error[i], " (prob.: ", x.probability[i],")")
    end
    println(io, "Integrand evaluations: ", x.neval)
    println(io, "Number of subregions:  ", x.nregions)
    print_fail(io, x)
end

@inline function dointegrate(x::T) where {T<:Integrand}
    # Choose the integrand function wrapper based on the value of `nvec`.  This function is
    # called only once, so the overhead of the following if should be negligible.
    if x.nvec == 1
        integrand = integrand_ptr(x.func)
    else
        integrand = integrand_ptr_nvec(x.func)
    end
    integral  = Vector{Cdouble}(undef, x.ncomp)
    error     = Vector{Cdouble}(undef, x.ncomp)
    prob      = Vector{Cdouble}(undef, x.ncomp)
    neval     = Ref{Int64}(0)
    fail      = Ref{Cint}(0)
    nregions  = Ref{Cint}(0)
    dointegrate!(x, integrand, integral, error, prob, neval, fail, nregions)
    return Integral{T}(integral, error, prob, neval[], fail[], nregions[])
end

### Other functions, not exported
function cores(n::Integer, p::Integer)
    ccall((:cubacores, libcuba), Ptr{Cvoid}, (Ptr{Cint}, Ptr{Cint}), Ref{Cint}(n), Ref{Cint}(p))
    return 0
end

function accel(n::Integer, p::Integer)
    ccall((:cubaaccel, libcuba), Ptr{Cvoid}, (Ptr{Cint}, Ptr{Cint}), Ref{Cint}(n), Ref{Cint}(p))
    return 0
end

function init(f::Ptr{Cvoid}, arg::Ptr{Cvoid})
    ccall((:cubainit, libcuba), Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}), f, arg)
    return 0
end

function exit(f::Ptr{Cvoid}, arg::Ptr{Cvoid})
    ccall((:cubaexit, libcuba), Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}), f, arg)
    return 0
end

end # module
