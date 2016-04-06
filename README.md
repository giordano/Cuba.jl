# Cuba

[![Build Status](https://travis-ci.org/giordano/Cuba.jl.svg?branch=master)](https://travis-ci.org/giordano/Cuba.jl) [![Coverage Status](https://coveralls.io/repos/github/giordano/Cuba.jl/badge.svg?branch=master)](https://coveralls.io/github/giordano/Cuba.jl?branch=master) [![codecov.io](https://codecov.io/github/giordano/Cuba.jl/coverage.svg?branch=master)](https://codecov.io/github/giordano/Cuba.jl?branch=master) [![Cuba](http://pkg.julialang.org/badges/Cuba_0.4.svg)](http://pkg.julialang.org/?pkg=Cuba) [![Cuba](http://pkg.julialang.org/badges/Cuba_0.5.svg)](http://pkg.julialang.org/?pkg=Cuba)

Introduction
------------

`Cuba.jl` is a library for multidimensional numerical integration with different
algorithms in [Julia](http://julialang.org/).

This is just a Julia wrapper around the C
[Cuba library](http://www.feynarts.de/cuba/) by **Thomas Hahn**.  All the
credits goes to him for the underlying functions, blame me for any problem with
the Julia interface.  Feel free to report bugs and make suggestions at
https://github.com/giordano/Cuba.jl/issues.

All algorithms provided by Cuba library are supported in `Cuba.jl`:
* `Vegas` (type: Monte Carlo; variance reduction with importance sampling)
* `Suave` (type: Monte Carlo; variance reduction with globally adaptive
  subdivision + importance sampling)
* `Divonne` (type: Monte Carlo or deterministic; variance reduction with
  stratified sampling, aided by methods from numerical optimization)
* `Cuhre` (type: deterministic; variance reduction with globally adaptive
  subdivision)

For more details on the algorithms see the manual included in Cuba library and
available at in `deps/cuba-shared-object/cuba.pdf` after successful installation
of `Cuba.jl`.

**Note:** This package works on the same operating systems for which the Cuba
library is available, i.e. GNU/Linux and OS X systems.

Installation
------------

`Cuba.jl` is available for Julia 0.4 and later versions, and can be installed
with
[Julia built-in package manager](http://docs.julialang.org/en/stable/manual/packages/).
In a Julia session run the command

```julia
julia> Pkg.add("Cuba")
```

You may need to update your package list with `Pkg.update()` in order to get the
latest version of `Cuba.jl`.


Usage
-----

After installing the package, run

``` julia
using Cuba
```

or put this command into your Julia script.

`Cuba.jl` provides 4 functions to integrate, one for each algorithm:

``` julia
Vegas(function, ndim, ncomp[, keywords...])
Suave(function, ndim, ncomp[, keywords...])
Divonne(function, ndim, ncomp[, keywords...])
Cuhre(function, ndim, ncomp[, keywords...])
```

Mandatory arguments are:

* `function`: the name of the function to be integrated
* `ndim`: the number of dimensions of the integral
* `ncomp`: the number of components of the integrand

The `function` must be of this type:

``` julia
function integrand(ndim::Cint, xx::Ptr{Cdouble}, ncomp::Cint, ff::Ptr{Cdouble},
                   userdata::Ptr{Void})
    # Take arrays from "xx" and "ff" pointers.
    x = pointer_to_array(xx, (ndim,))
    f = pointer_to_array(ff, (ncomp,))
	# Do calculations on "f" here
	#   ...
    # Store back the results to "ff"
    ff = pointer_from_objref(f)
return Cint(0)::Cint
end
```

Note that `xx` and `ff` arguments are passed as pointers, so you have to
translate them to Julia objects before actually performing calculations, and
finally store the results into `ff`.

All other arguments listed in Cuba documentation can be passed as optional
keywords.

The integrating functions `Vegas`, `Suave`, `Divonne`, and `Cuhre` return the
6-tuple

``` julia
(integral, error, probability, neval, fail, nregions)
```

The first three terms of the tuple are arrays with length `ncomp`, the last
three ones are scalars.  In particular, if you assign the output of integration
functions to the variable named `result`, you can access the value of the `i`-th
component of the integral with `result[1][i]` and the associated error with
`result[2][i]`.  The details of other quantities can be found in Cuba manual.

More extended documentation of `Cuba.jl` will come later.

**Note:** admittedly, this user interface is not REPL-friendly, help on
improving is welcome.

Example
-------

Here is an example of a 3-component integral in 3D space (so `ndim=3` and
`ncomp=3`) using the integrand function tested in `test/runtests.jl`:

``` julia
using Cuba

function func(ndim::Cint, xx::Ptr{Cdouble}, ncomp::Cint, ff::Ptr{Cdouble},
              userdata::Ptr{Void}=USERDATA_DEF)
    x = pointer_to_array(xx, (ndim,))
    f = pointer_to_array(ff, (ncomp,))
    f[1] = sin(x[1])*cos(x[2])*exp(x[3])
    f[2] = exp(-(x[1]^2 + x[2]^2 + x[3]^2))
    f[3] = 1/(1 - x[1]*x[2]*x[3])
    ff = pointer_from_objref(f)
    return Cint(0)::Cint
end

result = Cuhre(func, 3, 3, epsabs=1e-12, epsrel=1e-10)
println("Results of Cuba:")
println("Component 1: ", result[1][1], " ± ", result[2][1])
println("Component 2: ", result[1][2], " ± ", result[2][2])
println("Component 3: ", result[1][3], " ± ", result[2][3])
println("Exact results:")
println("Component 1: ", (e-1)*(1-cos(1))*sin(1))
println("Component 2: ", (sqrt(pi)*erf(1)/2)^3)
println("Component 3: ", zeta(3))
```

This is the output

```
Results of Cuba:
Component 1: 0.6646696797813739 ± 1.0050367631018485e-13
Component 2: 0.4165383858806454 ± 2.932866749838454e-11
Component 3: 1.2020569031649702 ± 1.1958522385908214e-10
Exact results:
Component 1: 0.6646696797813771
Component 2: 0.41653838588663805
Component 3: 1.2020569031595951
```

Performance
-----------

`Cuba.jl` has performances comparable with (if not slightly better than) an
equivalent native C code based on Cuba library when `CUBACORES` environment
variable is set to `0`.  This is the result of running the benchmark test
present in `test` directory.

```
$ CUBACORES=0 julia -f test/benchmark.jl
  [...]
INFO: Performance of Cuba.jl:
  0.338188 seconds (6.05 M allocations: 184.480 MB, 2.55% gc time)
  0.659029 seconds (6.00 M allocations: 183.107 MB, 1.21% gc time)
  0.384740 seconds (6.00 M allocations: 183.165 MB, 1.99% gc time)
  0.304015 seconds (6.00 M allocations: 183.129 MB, 2.86% gc time)
INFO: Performance of Cuba C Library:
  0.346084 seconds (Vegas)
  0.661870 seconds (Suave)
  0.378409 seconds (Divonne)
  0.306150 seconds (Cuhre)
```

Native C Cuba Library outperforms `Cuba.jl` when higher values of `CUBACORES` are
used, for example:

```
$ CUBACORES=1 julia -f test/benchmark.jl
  [...]
INFO: Performance of Cuba.jl:
  0.341448 seconds (6.05 M allocations: 184.480 MB, 2.60% gc time)
  0.660508 seconds (6.00 M allocations: 183.107 MB, 1.19% gc time)
  0.384731 seconds (6.00 M allocations: 183.165 MB, 2.01% gc time)
  0.302969 seconds (6.00 M allocations: 183.129 MB, 2.88% gc time)
INFO: Performance of Cuba C Library:
  0.119161 seconds (Vegas)
  0.608906 seconds (Suave)
  0.156459 seconds (Divonne)
  0.085269 seconds (Cuhre)
```

In `Cuba.jl`, the number of cores and accelerators to be used can be set with

``` julia
Cuba.cores(n, p)
Cuba.accel(n, p)
```

See Cuba manual for more information on these functions.

Related projects
----------------

Another Julia package for multidimenensional numerical integration is available:
[Cubature.jl](https://github.com/stevengj/Cubature.jl), by Steven G. Johnson.
Differently from `Cuba.jl`, this works on GNU/Linux, OS X and Windows as well.

License
-------

The Cuba.jl package is licensed under the GNU Lesser General Public License, the
same as [Cuba library](http://www.feynarts.de/cuba/).  The original author is
Mosè Giordano.  If you use this library for your work, please credit Thomas Hahn
(citable papers about Cuba library:
http://adsabs.harvard.edu/abs/2005CoPhC.168...78H and
http://adsabs.harvard.edu/abs/2015JPhCS.608a2066H).
