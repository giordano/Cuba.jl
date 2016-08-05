v0.1.1 (2016-08-05)
===================

Bug Fixes
------------

* Fix warnings in Julia 0.5.

v0.1.0 (2016-06-08)
===================

New Features
------------

* Module precompilation enabled

v0.0.5 (2016-04-12)
===================

New Features
------------

* User interface greatly simplified (thanks to Steven G. Johnson).  This change
  is **backward incompatible**.  See documentation for details

v0.0.4 (2016-04-10)
===================

New Features
------------

* New complete documentation, available at http://cubajl.readthedocs.org/ and
  locally in `docs/` directory

Breaking Changes
----------------

* `verbose` keyword renamed to `flags`

Bug Fixes
---------

* Number of cores fixed to 0 to avoid crashes when Julia has more than 1 process
* In `Cuhre` and `Divonne`, force `ndim` to be 2 when user sets it to 1

v0.0.3 (2016-04-06)
===================

New Features
------------

* Add `cores`, `accel`, `init`, `exit` function.  They will likely not be much
  useful for most users, so they are not exported nor documented.  See Cuba
  manual for information

Breaking Changes
----------------

* Make `ndim` and `ncomp` arguments mandatory

Bug Fixes
---------

* Fix build script

v0.0.2 (2016-04-04)
===================

Bug Fixes
---------

* Fix path of libcuba

v0.0.1 (2016-04-04)
===================

* First release
