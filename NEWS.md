History of Cuba.jl
==================

v2.0.0 (2019-02-27)
-------------------

### Breaking Changes

* Support for Julia 0.7 was dropped.

### New Features

* When using the package in the REPL, the result of the integration now has a
  more informative description about the error flag.

v1.0.0 (2018-08-17)
-------------------

### Breaking Changes

* Support for Julia 0.6 was dropped.
* Keyword arguments `reltol` and `abstol` are now called `rtol` and `atol`,
  respectively, to match keywords in `isapprox` function.
* Deprecated functions `llvegas`, `llsuave`, `lldivonne`, and `llcuhre` have
  been removed.

### New Features

* The build script has been updated, now the package supports Linux-musl and
  FreeBSD systems.

v0.5.0 (2018-05-15)
-------------------

### New Features

* The package now
  uses
  [`BinaryProvider.jl`](https://github.com/JuliaPackaging/BinaryProvider.jl) to
  install a pre-built version of the Cuba library on all platforms.

### Breaking Changes

* Support for Julia 0.5 was dropped
* The default value of argument `ndim` has been changed to 2 in `divonne` and
  `cuhre`.  These algorithms require the number of dimensions to be at least 2.
  Now setting `ndim` to 1 throws an error.  Your code will not be affected by
  this change if you did not explicitely set `ndim`.  See
  issue [#14](https://github.com/giordano/Cuba.jl/issues/14).

v0.4.0 (2017-07-08)
-------------------

### Breaking Changes

* Now `vegas`, `suave`, `divonne`, and `cuhre` wrap the 64-bit integers
  functions.  The 32-bit integers functions are no more available.  `llvegas`,
  `llsuave`, `lldivonne`, `llcuhre` are deprecated and will be removed at some
  point in the future.  This change reduces confusion about the function to use.

### New Features

* Now it is possible to vectorize a function in order to speed up its evaluation
  (see issue [#10](https://github.com/giordano/Cuba.jl/issues/10) and
  PR [#11](https://github.com/giordano/Cuba.jl/pull/11)).
* The result of integration is wrapped in an `Integral` object.  This is not a
  breaking change because its fields can be iterated over like a tuple, exactly
  as before.

v0.3.1 (2017-05-02)
-------------------

### Improvements

* Small performance improvements by avoiding dynamic dispatch in callback
  ([#6](https://github.com/giordano/Cuba.jl/pull/6)).  No user visible change.

v0.3.0 (2017-01-24)
-------------------

### Breaking Changes

* Support for Julia 0.4 was dropped
* Integrators functions with uppercase names were removed.  They were deprecated
  in v0.2.0

### New Features

* New 64-bit integers functions `llvegas`, `llsuave`, `lldivonne`, `llcuhre` are
  provided.  They should be used in cases where convergence is not reached
  within the ordinary 32-bit integer
  range ([#4](https://github.com/giordano/Cuba.jl/issues/4))

v0.2.0 (2016-10-15)
-------------------

This release faces some changes to the user interface.  Be aware of them when
upgrading.

### New Features

* `ndim` and `ncomp` arguments can be omitted.  In that case they default to 1.
  This change is not breaking, old syntax will continue to work.

### Breaking Changes

All integrator functions and some optional keywords have been renamed for more
consistency with the Julia environment.  Here is the detailed list:

* Integrators functions have been renamed to lowercase name: `Vegas` to `vegas`,
  `Suave` to `suave`, `Divonne` to `divonne`, `Cuhre` to `cuhre`.  The uppercase
  variants are still available but deprecated, they will be removed at some
  point in the future.
* Optional keywords changes: `epsabs` to `abstol`, `epsrel` to `reltol`,
  `maxeval` to `maxevals`, `mineval` to `minevals`.

v0.1.4 (2016-08-21)
-------------------

* A new version of Cuba library is downloaded to be compiled on GNU/Linux and
  Mac OS systems.  There have been only small changes for compatibility with
  recent GCC versions, no actual change to the library nor to the Julia wrapper.
  Nothing changed for Windows users.

v0.1.3 (2016-08-11)
-------------------

### New Features

* A tagged version of Cuba library is now downloaded when building the package.
  This ensures reproducibility of the results of a given `Cuba.jl` release.

v0.1.2 (2016-08-11)
-------------------

### New Features

* Windows (`i686` and `x86_64` architectures) supported
  ([#2](https://github.com/giordano/Cuba.jl/issues/2))

v0.1.1 (2016-08-05)
-------------------

### Bug Fixes

* Fix warnings in Julia 0.5

v0.1.0 (2016-06-08)
-------------------

### New Features

* Module precompilation enabled

v0.0.5 (2016-04-12)
-------------------

### New Features

* User interface greatly simplified (thanks to Steven G. Johnson;
  [#3](https://github.com/giordano/Cuba.jl/issues/3)).  This change is
  **backward incompatible**.  See documentation for details

v0.0.4 (2016-04-10)
-------------------

### New Features

* New complete documentation, available at http://cubajl.readthedocs.org/ and
  locally in `docs/` directory

### Breaking Changes

* `verbose` keyword renamed to `flags`

### Bug Fixes

* Number of cores fixed to 0 to avoid crashes when Julia has more than 1 process
* In `Cuhre` and `Divonne`, force `ndim` to be 2 when user sets it to 1

v0.0.3 (2016-04-06)
-------------------

### New Features

* Add `cores`, `accel`, `init`, `exit` function.  They will likely not be much
  useful for most users, so they are not exported nor documented.  See Cuba
  manual for information

### Breaking Changes

* Make `ndim` and `ncomp` arguments mandatory

### Bug Fixes

* Fix build script

v0.0.2 (2016-04-04)
-------------------

### Bug Fixes

* Fix path of libcuba

v0.0.1 (2016-04-04)
-------------------

* First release
