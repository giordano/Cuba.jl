# Cuba

[![Build Status](https://travis-ci.org/giordano/Cuba.jl.svg?branch=master)](https://travis-ci.org/giordano/Cuba.jl)

Library for multidimensional numerical integration with different algorithms.
Sampling of the integrand is multi-threaded.

This is just a Julia wrapper around the C
[Cuba library](http://www.feynarts.de/cuba/) by **Thomas Hahn**.  All the credits
goes to him for the underlying functions, blame me for any problem with the
Julia interface.  Feel free to report bugs at
https://github.com/giordano/Cuba.jl/issues.

Note: Cuba library works only on GNU/Linux and OS X systems.  Currently, this
Julia package works only on GNU/Linux, support for OS X will come soon-ish.

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
