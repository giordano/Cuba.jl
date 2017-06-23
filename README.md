# Cuba.jl

| **Documentation**                       | [**Package Evaluator**][pkgeval-link] | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![][pkg-0.5-img]][pkg-0.5-url]       | [![Build Status][travis-img]][travis-url] | [![][coveral-img]][coveral-url] |
| [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.6-img]][pkg-0.6-url]       | [![Build Status][appvey-img]][appvey-url] | [![][codecov-img]][codecov-url] |

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

For more details on the algorithms see the manual included in Cuba library and
available in `deps/cuba-julia/cuba.pdf` after successful installation
of `Cuba.jl`.

Integration is performed on the n-dimensional unit hypercube $[0, 1]^n$.  If you
want to compute an integral over a different set, you have to scale the
integrand function in order to have an equivalent integral on $[0, 1]^n$.  For
example, recall that in one dimension

```
∫_a^b dx f[x] → ∫_0^1 dy f[a + (b - a) y] (b - a)
```

where the final `(b - a)` is the one-dimensional version of the Jacobian.  This
generalizes straightforwardly to more than one dimension.

`Cuba.jl` is available for GNU/Linux, Mac OS, and Windows (`i686` and `x86_64`
architectures).

Installation
------------

`Cuba.jl` is available for Julia 0.5 and later versions, and can be installed
with
[Julia built-in package manager](http://docs.julialang.org/en/stable/manual/packages/).
In a Julia session run the commands

```julia
julia> Pkg.update()
julia> Pkg.add("Cuba")
```

Installation script on GNU/Linux and Mac OS systems will download Cuba Library
source code and build the Cuba shared object.  In order to accomplish this task
a C compiler is needed.  Instead, on Windows a prebuilt version of the library
is downloaded.

Older versions are also available for Julia 0.4.

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

* `ndim`: the number of dimensions of the integration domain.  Defaults to 1
* `ncomp`: the number of components of the integrand.  Defaults to 1

`ndim` and `ncomp` arguments must appear in this order, so you cannot omit
`ndim` but not `ncomp`.  `integrand` should be a function `integrand(x, f)`
taking two arguments:

- the input vector `x` of length `ndim`
- the output vector `f` of length `ncomp`, used to set the value of each
  component of the integrand at point `x`

Also
[anonymous functions](http://docs.julialang.org/en/stable/manual/functions/#anonymous-functions)
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

**Note:** if you used `Cuba.jl` until version 0.0.4, be aware that the user
interface has been reworked in version 0.0.5 in a backward incompatible way.

### Documentation ###

A more detailed manual of `Cuba.jl`, with many complete examples, is available
at http://cubajl.readthedocs.io/.  You can also download the latest PDF version
from https://media.readthedocs.org/pdf/cubajl/latest/cubajl.pdf.

Example
-------

Here is an example of a 3-component integral in 3D space (so `ndim=3` and
`ncomp=3`) using the integrand function tested in `test/runtests.jl`:

``` julia
using Cuba

function func(x, f)
    f[1] = sin(x[1])*cos(x[2])*exp(x[3])
    f[2] = exp(-(x[1]^2 + x[2]^2 + x[3]^2))
    f[3] = 1/(1 - x[1]*x[2]*x[3])
end

result = cuhre(func, 3, 3, abstol=1e-12, reltol=1e-10)
println("Results of Cuba:")
for i=1:3; println("Component $i: ", result[1][i], " ± ", result[2][i]); end
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

`Cuba.jl` cannot ([yet?](https://github.com/giordano/Cuba.jl/issues/1)) take
advantage of parallelization capabilities of Cuba Library.  Nonetheless, it has
performances competitive with equivalent native C or Fortran codes based on Cuba
library when `CUBACORES` environment variable is set to `0` (i.e.,
multithreading is disabled).  The following is the result of running the
benchmark present in `test` directory on a 64-bit GNU/Linux system running Julia
0.7.0-DEV.363 (commit 6071f1a02e) equipped with an Intel(R) Core(TM) i7-4700MQ
CPU.  The C and FORTRAN 77 benchmark codes have been compiled with GCC 6.3.0.

```
$ CUBACORES=0 julia -e 'cd(Pkg.dir("Cuba")); include("test/benchmark.jl")'
INFO: Performance of Cuba.jl:
  0.271304 seconds (Vegas)
  0.579783 seconds (Suave)
  0.329504 seconds (Divonne)
  0.238852 seconds (Cuhre)
INFO: Performance of Cuba Library in C:
  0.319799 seconds (Vegas)
  0.619774 seconds (Suave)
  0.340317 seconds (Divonne)
  0.266906 seconds (Cuhre)
INFO: Performance of Cuba Library in Fortran:
  0.272000 seconds (Vegas)
  0.584000 seconds (Suave)
  0.308000 seconds (Divonne)
  0.232000 seconds (Cuhre)
```

Of course, native C and Fortran codes making use of Cuba Library outperform
`Cuba.jl` when higher values of `CUBACORES` are used, for example:

```
$ CUBACORES=1 julia -e 'cd(Pkg.dir("Cuba")); include("test/benchmark.jl")'
INFO: Performance of Cuba.jl:
  0.279524 seconds (Vegas)
  0.581078 seconds (Suave)
  0.327319 seconds (Divonne)
  0.241211 seconds (Cuhre)
INFO: Performance of Cuba Library in C:
  0.115113 seconds (Vegas)
  0.596503 seconds (Suave)
  0.152511 seconds (Divonne)
  0.085805 seconds (Cuhre)
INFO: Performance of Cuba Library in Fortran:
  0.108000 seconds (Vegas)
  0.604000 seconds (Suave)
  0.160000 seconds (Divonne)
  0.092000 seconds (Cuhre)
```

`Cuba.jl` internally fixes `CUBACORES` to 0 in order to prevent from forking
`julia` processes that would only slow down calculations eating up the memory,
without actually taking advantage of concurrency.  Furthemore, without this
measure, adding more Julia processes with `addprocs()` would only make the
program segfault.

Related projects
----------------

Another Julia package for multidimenensional numerical integration is available:
[Cubature.jl](https://github.com/stevengj/Cubature.jl), by Steven G. Johnson.

License
-------

The Cuba.jl package is licensed under the GNU Lesser General Public License, the
same as [Cuba library](http://www.feynarts.de/cuba/).  The original author is
Mosè Giordano.  If you use this library for your work, please credit Thomas Hahn
(citable papers about Cuba library:
http://adsabs.harvard.edu/abs/2005CoPhC.168...78H and
http://adsabs.harvard.edu/abs/2015JPhCS.608a2066H).



[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://cubajl.readthedocs.io/en/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://cubajl.readthedocs.io/en/stable/

[pkgeval-link]: http://pkg.julialang.org/?pkg=Cuba

[pkg-0.5-img]: http://pkg.julialang.org/badges/Cuba_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/detail/Cuba.html
[pkg-0.6-img]: http://pkg.julialang.org/badges/Cuba_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/detail/Cuba.html

[travis-img]: https://travis-ci.org/giordano/Cuba.jl.svg?branch=master
[travis-url]: https://travis-ci.org/giordano/Cuba.jl

[appvey-img]: https://ci.appveyor.com/api/projects/status/ivqy72upfjxplbcn/branch/master?svg=true
[appvey-url]: https://ci.appveyor.com/project/giordano/cuba-jl

[coveral-img]: https://coveralls.io/repos/github/giordano/Cuba.jl/badge.svg?branch=master
[coveral-url]: https://coveralls.io/github/giordano/Cuba.jl?branch=master

[codecov-img]: https://codecov.io/gh/giordano/Cuba.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/giordano/Cuba.jl
