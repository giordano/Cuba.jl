### suave.jl --- Integrate with Suave method

struct Suave{T, D} <: Integrand{T}
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
    nnew::Int64
    nmin::Int64
    flatness::Cdouble
    statefile::String
    spin::Ptr{Cvoid}
end

@inline function dointegrate!(x::Suave{T, D}, integrand, integral,
                              error, prob, neval, fail, nregions) where {T, D}

    userdata = ismissing(x.userdata) ? x.func : (x.func, x.userdata)

    ccall((:llSuave, libcuba), Cdouble,
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
          x.ndim, x.ncomp, integrand, userdata, x.nvec,
          x.rtol, x.atol, x.flags, x.seed, x.minevals, x.maxevals,
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

* `userdata`
* `nvec`
* `rtol`
* `atol`
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
               nvec::Integer=NVEC, rtol::Real=RTOL, atol::Real=ATOL,
               flags::Integer=FLAGS, seed::Integer=SEED,
               minevals::Real=MINEVALS, maxevals::Real=MAXEVALS,
               nnew::Integer=NNEW, nmin::Integer=NMIN, flatness::Real=FLATNESS,
               statefile::AbstractString=STATEFILE, spin::Ptr{Cvoid}=SPIN,
               reltol=missing, abstol=missing, userdata=missing) where {T}
    atol_,rtol_ = tols(atol,rtol,abstol,reltol)
    # See <https://github.com/giordano/Cuba.jl/issues/27>.
    max_maxevals = typemax(Int64) ÷ 2
    maxevals > max_maxevals && throw(ArgumentError("`maxevals` can't be larger than $(max_maxevals) in `suave`, found $(maxevals)"))
    return dointegrate(Suave(integrand, userdata, ndim, ncomp, Int64(nvec), Cdouble(rtol_),
                             Cdouble(atol_), flags, seed, trunc(Int64, minevals),
                             trunc(Int64, maxevals), Int64(nnew), Int64(nmin),
                             Cdouble(flatness), String(statefile), spin))
end
