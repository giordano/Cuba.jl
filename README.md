# Cuba.jl

| **Documentation**                       | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][travis-img]][travis-url] | [![][coveral-img]][coveral-url] |
| [![][docs-latest-img]][docs-latest-url] | [![Build Status][appvey-img]][appvey-url] | [![][codecov-img]][codecov-url] |

Introduction
------------

`Cuba.jl` is a library for multidimensional numerical integration with different
algorithms in [Julia](http://julialang.org/).

This is just a Julia wrapper around the C
[Cuba library](http://www.feynarts.de/cuba/), version 4.2, by **Thomas Hahn**.
All the credits goes to him for the underlying functions, blame me for any
problem with the Julia interface.  Feel free to report bugs and make suggestions
at https://github.com/giordano/Cuba.jl/issues.

All algorithms provided by Cuba library are supported in `Cuba.jl`:

* `vegas` (type: Monte Carlo; variance reduction with importance sampling)
* `suave` (type: Monte Carlo; variance reduction with globally adaptive
  subdivision + importance sampling)
* `divonne` (type: Monte Carlo or deterministic; variance reduction with
  stratified sampling, aided by methods from numerical optimization)
* `cuhre` (type: deterministic; variance reduction with globally adaptive
  subdivision)

Integration is performed on the n-dimensional unit hypercube [0, 1]^n.  For more
details on the algorithms see the manual included in Cuba library and available
in `deps/usr/share/cuba.pdf` after successful installation of `Cuba.jl`.

`Cuba.jl` is available for GNU/Linux, FreeBSD, Mac OS, and Windows (`i686` and
`x86_64` architectures).

Installation
------------

The latest version of `Cuba.jl` is available for Julia 1.0 and later versions,
and can be installed with [Julia built-in package
manager](https://docs.julialang.org/en/v1/stdlib/Pkg/). In a Julia session run
the commands

```julia
pkg> update
pkg> add Cuba
```

Older versions are also available for Julia 0.4-0.7.

Usage
-----

After installing the package, run

``` julia
using Cuba
```

or put this command into your Julia script.

`Cuba.jl` provides the following functions to integrate:

``` julia
vegas(integrand, ndim, ncomp[; keywords...])
suave(integrand, ndim, ncomp[; keywords...])
divonne(integrand, ndim, ncomp[; keywords...])
cuhre(integrand, ndim, ncomp[; keywords...])
```

These functions wrap the 64-bit integers functions provided by the Cuba library.

The only mandatory argument is:

* `function`: the name of the function to be integrated

Optional positional arguments are:

* `ndim`: the number of dimensions of the integration domain.  Defaults to 1 in
  `vegas` and `suave`, to 2 in `divonne` and `cuhre`.  Note: `ndim` must be
  at least 2 with the latest two methods.
* `ncomp`: the number of components of the integrand.  Defaults to 1

`ndim` and `ncomp` arguments must appear in this order, so you cannot omit
`ndim` but not `ncomp`.  `integrand` should be a function `integrand(x, f)`
taking two arguments:

- the input vector `x` of length `ndim`
- the output vector `f` of length `ncomp`, used to set the value of each
  component of the integrand at point `x`

Also
[anonymous functions](https://docs.julialang.org/en/v1/manual/functions/#man-anonymous-functions-1)
are allowed as `integrand`.  For those familiar with
[`Cubature.jl`](https://github.com/stevengj/Cubature.jl) package, this is the
same syntax used for integrating vector-valued functions.

For example, the integral

```
∫_0^1 cos(x) dx = sin(1) = 0.8414709848078965
```

can be computed with one of the following commands

``` julia
julia> vegas((x, f) -> f[1] = cos(x[1]))
Component:
 1: 0.8414910005259612 ± 5.2708169787342156e-5 (prob.: 0.028607201258072673)
Integrand evaluations: 13500
Fail:                  0
Number of subregions:  0

julia> suave((x, f) -> f[1] = cos(x[1]))
Component:
 1: 0.84115236906584 ± 8.357995609919512e-5 (prob.: 1.0)
Integrand evaluations: 22000
Fail:                  0
Number of subregions:  22

julia> divonne((x, f) -> f[1] = cos(x[1]))
Component:
 1: 0.841468071955942 ± 5.3955070531551656e-5 (prob.: 0.0)
Integrand evaluations: 1686
Fail:                  0
Number of subregions:  14

julia> cuhre((x, f) -> f[1] = cos(x[1]))
Component:
 1: 0.8414709848078966 ± 2.2204460420128823e-16 (prob.: 3.443539937576958e-5)
Integrand evaluations: 195
Fail:                  0
Number of subregions:  2
```

The integrating functions `vegas`, `suave`, `divonne`, and `cuhre` return an
`Integral` object whose fields are

``` julia
integral :: Vector{Float64}
error    :: Vector{Float64}
probl    :: Vector{Float64}
neval    :: Int64
fail     :: Int32
nregions :: Int32
```

The first three fields are vectors with length `ncomp`, the last three ones are
scalars.  The `Integral` object can also be iterated over like a tuple.  In
particular, if you assign the output of integration functions to the variable
named `result`, you can access the value of the `i`-th component of the integral
with `result[1][i]` or `result.integral[i]` and the associated error with
`result[2][i]` or `result.error[i]`.  The details of other quantities can be
found in Cuba manual.

All other arguments listed in Cuba documentation can be passed as optional
keywords.

### Documentation ###

A more detailed manual of `Cuba.jl`, with many complete examples, is available
at https://giordano.github.io/Cuba.jl/stable/.

Related projects
----------------

There are other Julia packages for multidimenensional numerical integration:

* [`Cubature.jl`](https://github.com/stevengj/Cubature.jl)
* [`HCubature.jl`](https://github.com/stevengj/HCubature.jl)
* [`NIntegration.jl`](https://github.com/pabloferz/NIntegration.jl)

License
-------

The Cuba.jl package is licensed under the GNU Lesser General Public License, the
same as [Cuba library](http://www.feynarts.de/cuba/).  The original author is
Mosè Giordano.  If you use this library for your work, please credit Thomas Hahn
(citable papers about Cuba library:
http://adsabs.harvard.edu/abs/2005CoPhC.168...78H and
http://adsabs.harvard.edu/abs/2015JPhCS.608a2066H).



[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://giordano.github.io/Cuba.jl/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://giordano.github.io/Cuba.jl/stable/

[travis-img]: https://travis-ci.org/giordano/Cuba.jl.svg?branch=master
[travis-url]: https://travis-ci.org/giordano/Cuba.jl

[appvey-img]: https://ci.appveyor.com/api/projects/status/ivqy72upfjxplbcn/branch/master?svg=true
[appvey-url]: https://ci.appveyor.com/project/giordano/cuba-jl

[coveral-img]: https://coveralls.io/repos/github/giordano/Cuba.jl/badge.svg?branch=master
[coveral-url]: https://coveralls.io/github/giordano/Cuba.jl?branch=master

[codecov-img]: https://codecov.io/gh/giordano/Cuba.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/giordano/Cuba.jl
