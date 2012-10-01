Major:

* Implement the documentation component. (I don't grok what's going on in
  this makefiles. There is a rudimentary doxygen call, though.)

* Add the rest of the Bug tests

* Test the fparser 4.5 functionality.

* The following features have to be written:
    * MUMPS
    * SLEPC
    * HDF5
    * TECPLOT

* Write a wrapper around cmake. Something like
    * download and 'install' cmake if necessary
    $ make config
    $ make
    ...

* Update the online documentation on

    - How to install deal.II

    - How to interface with external libraries (Unify this...)

* Port the testsuite to cmake

* Build system for parameter_gui.

* Finish install targets:

   - share files (documentation/manpage/etc.)
   - same for mesh_conversion and parameter_gui

* Rename all DEAL_II_USE_* and HAVE_LIB* that refer to a feature to
  DEAL_II_WITH_*.


Minor:

* Get rid of DEAL_II_PATH in config.h and in the sources.

* HAVE_JN_ ?

* There is still a perl script for generating lapack_templates.h. Is this
  necessary.

* TODO: Update the p4est-setup.sh script...

