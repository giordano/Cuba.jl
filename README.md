# Cuba

[![Build Status](https://travis-ci.org/giordano/Cuba.jl.svg?branch=master)](https://travis-ci.org/giordano/Cuba.jl) [![Coverage Status](https://coveralls.io/repos/github/giordano/Cuba.jl/badge.svg?branch=master)](https://coveralls.io/github/giordano/Cuba.jl?branch=master) [![codecov.io](https://codecov.io/github/giordano/Cuba.jl/coverage.svg?branch=master)](https://codecov.io/github/giordano/Cuba.jl?branch=master) [![Cuba](http://pkg.julialang.org/badges/Cuba_0.4.svg)](http://pkg.julialang.org/?pkg=Cuba) [![Cuba](http://pkg.julialang.org/badges/Cuba_0.5.svg)](http://pkg.julialang.org/?pkg=Cuba)

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

* `Vegas` (type: Monte Carlo; variance reduction with importance sampling)
* `Suave` (type: Monte Carlo; variance reduction with globally adaptive
  subdivision + importance sampling)
* `Divonne` (type: Monte Carlo or deterministic; variance reduction with
  stratified sampling, aided by methods from numerical optimization)
* `Cuhre` (type: deterministic; variance reduction with globally adaptive
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

**Note:** This package has been tested only on GNU/Linux and OS X systems.
Trying to install on Windows will likely fail, please report at
https://github.com/giordano/Cuba.jl/issues/2 if you manage to install on this
system.

Installation
------------

`Cuba.jl` is available for Julia 0.4 and later versions, and can be installed
with
[Julia built-in package manager](http://docs.julialang.org/en/stable/manual/packages/).
In a Julia session run the command

```julia
julia> Pkg.add("Cuba")
```

Installation script will download Cuba Library source code and build the Cuba
shared object.  In order to accomplish this task a C compiler is needed.

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
Vegas(integrand, ndim, ncomp[, keywords...])
Suave(integrand, ndim, ncomp[, keywords...])
Divonne(integrand, ndim, ncomp[, keywords...])
Cuhre(integrand, ndim, ncomp[, keywords...])
```

Mandatory arguments are:

* `function`: the name of the function to be integrated
* `ndim`: the number of dimensions of the integral
* `ncomp`: the number of components of the integrand

`integrand` should be a function `integrand(x, f)` taking two arguments:

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
∫_0^1 cos(x) dx
```

can be computed with one of the following lines

``` julia
Vegas((x,f)->f[1]=cos(x[1]), 1, 1)
Suave((x,f)->f[1]=cos(x[1]), 1, 1)
Divonne((x,f)->f[1]=cos(x[1]), 1, 1)
Cuhre((x,f)->f[1]=cos(x[1]), 1, 1)
```

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

All other arguments listed in Cuba documentation can be passed as optional
keywords.  More extended documentation of `Cuba.jl` is available at
http://cubajl.readthedocs.org/.  You can also download the PDF version of the
manual from https://media.readthedocs.org/pdf/cubajl/latest/cubajl.pdf.

**Note:** if you used `Cuba.jl` until version 0.4, be aware that the user
interface has been reworked in version 0.5 in a backward incompatible way.

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

result = Cuhre(func, 3, 3, epsabs=1e-12, epsrel=1e-10)
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

`Cuba.jl` cannot (yet?) take advantage of parallelization capabilities of Cuba
Library.  Nonetheless, it has performances comparable with equivalent native C
or Fortran codes based on Cuba library when `CUBACORES` environment variable is
set to `0` (i.e., multithreading is disabled).  The following is the result of
running the benchmark present in `test` directory on a 64-bit GNU/Linux system
running Julia 0.4.3.  The C and FORTRAN 77 benchmark codes have been compiled
with GCC 5.3.1.

```
$ CUBACORES=0 julia -e 'cd(Pkg.dir("Cuba")); include("test/benchmark.jl")'
INFO: Performance of Cuba.jl:
  0.340635 seconds (Vegas)
  0.660305 seconds (Suave)
  0.391721 seconds (Divonne)
  0.305756 seconds (Cuhre)
INFO: Performance of Cuba Library in C:
  0.352429 seconds (Vegas)
  0.668258 seconds (Suave)
  0.380006 seconds (Divonne)
  0.305772 seconds (Cuhre)
INFO: Performance of Cuba Library in Fortran:
  0.328000 seconds (Vegas)
  0.660000 seconds (Suave)
  0.364000 seconds (Divonne)
  0.296000 seconds (Cuhre)

```

Of course, native C and Fortran codes making use of Cuba Library outperform
`Cuba.jl` when higher values of `CUBACORES` are used, for example:

```
$ CUBACORES=1 julia -e 'cd(Pkg.dir("Cuba")); include("test/benchmark.jl")'
INFO: Performance of Cuba.jl:
  0.342575 seconds (Vegas)
  0.660071 seconds (Suave)
  0.393213 seconds (Divonne)
  0.304569 seconds (Cuhre)
INFO: Performance of Cuba Library in C:
  0.118911 seconds (Vegas)
  0.614480 seconds (Suave)
  0.153015 seconds (Divonne)
  0.086997 seconds (Cuhre)
INFO: Performance of Cuba Library in Fortran:
  0.108000 seconds (Vegas)
  0.628000 seconds (Suave)
  0.144000 seconds (Divonne)
  0.084000 seconds (Cuhre)
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
Differently from `Cuba.jl`, this works on GNU/Linux, OS X and Windows as well.

License
-------

The Cuba.jl package is licensed under the GNU Lesser General Public License, the
same as [Cuba library](http://www.feynarts.de/cuba/).  The original author is
Mosè Giordano.  If you use this library for your work, please credit Thomas Hahn
(citable papers about Cuba library:
http://adsabs.harvard.edu/abs/2005CoPhC.168...78H and
http://adsabs.harvard.edu/abs/2015JPhCS.608a2066H).
