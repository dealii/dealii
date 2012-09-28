Major:

* Implement the documentation component. (I don't grok what's going on in
  this makefiles. There is a rudimentary doxygen call, though.)

* Add the rest of the Bug tests

* The following features have to be written:
    * MUMPS
    * SLEPC
    * HDF5
    * HSL
    * TECPLOT

* Write a wrapper around cmake. Something like
    $ make config
    $ make
    ...

* Update the online documentation on

    - How to install deal.II

    - How to interface with external libraries (Unify this...)

* Port the testsuite to cmake


Minor:

* There is still a perl script for generating lapack_templates.h. Is this
  necessary.

* What to do with ./lib/meshes ? I have moved it to ./contrib/meshes for
  the time being...

* The original decision to have a dedicated p4est/FAST and p4est/DEBUG
  structure is a bit unfortunate and therefor replaced by a FindP4EST
  module that searches libp4est.so and treats this as FAST and
  libp4est_debug.so or libp4est.g.so and treats this as debug library.
  TODO: Update the p4est-setup.sh script...

* Reenable support for metis versions 4.x ?!

* The contrib libraries are included too early. User CXX_FLAGS won't be
  included in the contrib object targets...
