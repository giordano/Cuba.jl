### divonne.jl --- Integrate with Divonne method

struct Divonne{T} <: Integrand{T}
    func::T
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

@inline function dointegrate!(x::Divonne{T}, integrand, integral,
                              error, prob, neval, fail, nregions) where {T}
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
          x.ndim, x.ncomp, integrand, x.func, x.nvec, x.rtol,
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
                 spin::Ptr{Cvoid}=SPIN, reltol=missing, abstol=missing) where {T}
    atol_,rtol_ = tols(atol,rtol,abstol,reltol)
    # Divonne requires "ndim" to be at least 2, even for an integral over a one
    # dimensional domain.  Instead, we don't prevent users from setting wrong
    # "ndim" values like 0 or negative ones.
    ndim >=2 || throw(ArgumentError("In Divonne ndim must be ≥ 2"))
    return dointegrate(Divonne(integrand, ndim, ncomp, Int64(nvec), Cdouble(rtol_),
                               Cdouble(atol_), flags, seed, trunc(Int64, minevals),
                               trunc(Int64, maxevals), key1, key2, key3, maxpass,
                               Cdouble(border), Cdouble(maxchisq), Cdouble(mindeviation),
                               Int64(ngiven), ldxgiven, xgiven, Int64(nextra),
                               peakfinder, String(statefile), spin))
end
