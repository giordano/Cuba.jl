using Cuba
using Base.Test

function func(ndim::Cint,
              xx::Ptr{Cdouble},
              ncomp::Cint,
              ff::Ptr{Cdouble},
              userdata::Ptr{Void}=USERDATA_DEF
              )
    x = pointer_to_array(xx, (ndim,))
    f = pointer_to_array(ff, (ncomp,))

    # Define components of the integrand function.
    f[1] = sin(x[1])*cos(x[2])*exp(x[3])
    f[2] = exp(-(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]))
    f[3] = 1/(1 - x[1]*x[2]*x[3])

    xx = pointer_from_objref(x)
    ff = pointer_from_objref(f)
    return Cint(0)::Cint
end

# Test results and make sure the estimation of error is exact.
let
    local result
    @time result = Vegas(func, ndim=3, ncomp=3, epsabs=1e-6)
    @test_approx_eq_eps result[1][1]   (e-1)*(1-cos(1))*sin(1)   result[2][1]
    @test_approx_eq_eps result[1][2]   (sqrt(pi)*erf(1)/2)^3     result[2][2]
    @test_approx_eq_eps result[1][3]   zeta(3)                   result[2][3]
end
