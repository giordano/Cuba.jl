### cuhre.jl --- Integrate with Cuhre method

struct Cuhre{T, D} <: Integrand{T}
    func::T
    userdata::D
    ndim::Int
    ncomp::Int
    nvec::Int64
    rtol::Cdouble
    atol::Cdouble
    flags::Int
    minevals::Int64
    maxevals::Int64
    key::Int
    statefile::String
    spin::Ptr{Cvoid}
end

@inline function dointegrate!(x::Cuhre{T, D}, integrand, integral,
                              error, prob, neval, fail, nregions) where {T, D}

    userdata = ismissing(x.userdata) ? x.func : (x.func, x.userdata)

    ccall((:llCuhre, libcuba), Cdouble,
          (Cint, # ndim
           Cint, # ncomp
           Ptr{Cvoid}, # integrand
           Any, # userdata
           Int64, # nvec
           Cdouble, # rtol
           Cdouble, # atol
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
          x.ndim, x.ncomp, integrand, userdata, x.nvec, x.rtol, x.atol,
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

* `userdata`
* `nvec`
* `rtol`
* `atol`
* `flags`
* `minevals`
* `maxevals`
* `key`
* `statefile`
* `spin`
"""
function cuhre(integrand::T, ndim::Integer=2, ncomp::Integer=1;
               nvec::Integer=NVEC, rtol::Real=RTOL, atol::Real=ATOL,
               flags::Integer=FLAGS, minevals::Real=MINEVALS,
               maxevals::Real=MAXEVALS, key::Integer=KEY,
               statefile::AbstractString=STATEFILE, spin::Ptr{Cvoid}=SPIN,
               abstol=missing, reltol=missing, userdata=missing) where {T}
    atol_,rtol_ = tols(atol,rtol,abstol,reltol)
    # Cuhre requires "ndim" to be at least 2, even for an integral over a one
    # dimensional domain.  Instead, we don't prevent users from setting wrong
    # "ndim" values like 0 or negative ones.
    ndim >=2 || throw(ArgumentError("In Cuhre ndim must be ≥ 2"))
    return dointegrate(Cuhre(integrand, userdata, ndim, ncomp, Int64(nvec), Cdouble(rtol_),
                             Cdouble(atol_), flags, trunc(Int64, minevals),
                             trunc(Int64, maxevals), key, String(statefile), spin))
end
