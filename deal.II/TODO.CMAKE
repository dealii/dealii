Major:

* Implement the documentation component. (I don't grok what's going on in
  this makefiles. There is a rudimentary doxygen call, though.)

* Add the rest of the Bug tests

* Test the fparser 4.5 functionality.

* The following features have to be written:
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

* Rename all DEAL_II_USE_* and HAVE_LIB* that refer to a feature to
  DEAL_II_WITH_*.

* Refactor lapack_templates.pl


Minor:

* Fixup the TODOs in
    source/lac/sparse_dircect.cc
    source/base/path_find.cc
    source/base/function_parser.cc

* HAVE_JN_ ?

* TODO: Update the p4est-setup.sh script...

