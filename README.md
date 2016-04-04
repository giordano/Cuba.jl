# Cuba

[![Build Status](https://travis-ci.org/giordano/Cuba.jl.svg?branch=master)](https://travis-ci.org/giordano/Cuba.jl) [![Coverage Status](https://coveralls.io/repos/github/giordano/Cuba.jl/badge.svg?branch=master)](https://coveralls.io/github/giordano/Cuba.jl?branch=master) [![codecov.io](https://codecov.io/github/giordano/Cuba.jl/coverage.svg?branch=master)](https://codecov.io/github/giordano/Cuba.jl?branch=master) [![Cuba](http://pkg.julialang.org/badges/Cuba_0.4.svg)](http://pkg.julialang.org/?pkg=Cuba) [![Cuba](http://pkg.julialang.org/badges/Cuba_0.5.svg)](http://pkg.julialang.org/?pkg=Cuba)

Introduction
------------

`Cuba.jl` is a library for multidimensional numerical integration with different
algorithms in [Julia](http://julialang.org/).  Sampling of the integrand is
multi-threaded.

This is just a Julia wrapper around the C
[Cuba library](http://www.feynarts.de/cuba/) by **Thomas Hahn**.  All the
credits goes to him for the underlying functions, blame me for any problem with
the Julia interface.  Feel free to report bugs at
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

Usage
-----

After installing the package, run

``` julia
using Cuba
```

or put this command into your Julia script.

Sorry, no documentation so far, it will be added when the user interface will
reach a final state.  For the time being have a look at `test/runtests.jl` to
see how you can solve integrals.

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
