### vegas.jl --- Integrate with Vegas method

struct Vegas{T, D} <: Integrand{T}
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
    nstart::Int64
    nincrease::Int64
    nbatch::Int64
    gridno::Int
    statefile::String
    spin::Ptr{Cvoid}
end

@inline function dointegrate!(x::Vegas{T, D}, integrand, integral,
                              error, prob, neval, fail, nregions) where {T, D}

    userdata = ismissing(x.userdata) ? x.func : (x.func, x.userdata)

    ccall((:llVegas, libcuba), Cdouble,
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
           Int64, # nstart
           Int64, # nincrease
           Int64, # nbatch
           Cint, # gridno
           Ptr{Cchar}, # statefile
           Ptr{Cvoid}, # spin
           Ptr{Int64}, # neval
           Ptr{Cint}, # fail
           Ptr{Cdouble}, # integral
           Ptr{Cdouble}, # error
           Ptr{Cdouble}),# prob
          # Input
          x.ndim, x.ncomp, integrand, userdata, x.nvec,
          x.rtol, x.atol, x.flags, x.seed, x.minevals, x.maxevals,
          x.nstart, x.nincrease, x.nbatch, x.gridno, x.statefile, x.spin,
          # Output
          neval, fail, integral, error, prob)
end

"""
    vegas(integrand, ndim=1, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions

Calculate integral of `integrand` over the unit hypercube in `ndim` dimensions
using Vegas algorithm.  `integrand` is a vectorial function with `ncomp`
components.  `ndim` and `ncomp` default to 1.

Accepted keywords:

* `userdata`
* `nvec`
* `rtol`
* `atol`
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
function vegas(integrand::T, ndim::Integer=1, ncomp::Integer=1;
               nvec::Integer=NVEC, rtol::Real=RTOL, atol::Real=ATOL,
               flags::Integer=FLAGS, seed::Integer=SEED,
               minevals::Real=MINEVALS, maxevals::Real=MAXEVALS,
               nstart::Integer=NSTART, nincrease::Integer=NINCREASE,
               nbatch::Integer=NBATCH, gridno::Integer=GRIDNO,
               statefile::AbstractString=STATEFILE, spin::Ptr{Cvoid}=SPIN,
               reltol=missing, abstol=missing, userdata=missing) where {T}
    atol_,rtol_ = tols(atol,rtol,abstol,reltol)
    return dointegrate(Vegas(integrand, userdata, ndim, ncomp, Int64(nvec), Cdouble(rtol_),
                             Cdouble(atol_), flags, seed, trunc(Int64, minevals),
                             trunc(Int64, maxevals), Int64(nstart),
                             Int64(nincrease), Int64(nbatch), gridno,
                             String(statefile), spin))
end
