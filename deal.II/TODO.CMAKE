Major:

* MPI dependencies?
  arpack, hdf5, petsc, p4est, trilinos, metis

* Add the rest of the Bug tests

* Test the fparser 4.5 functionality.

* The following features have to be written:
    * SLEPC
    * TECPLOT

* Write a wrapper around cmake. Something like
    $ make config
    $ make
    ...

* Update the online documentation on

    - How to install deal.II

    - How to interface with external libraries (Unify this...)

* Port the testsuite to cmake

* Rename all DEAL_II_USE_* and HAVE_LIB* that refer to a feature to
  DEAL_II_WITH_*.


Minor:

* (Maybe) re-enable validate-xrefs

* Fixup the TODOs in
    source/lac/sparse_dircect.cc
    source/base/path_find.cc
    source/base/function_parser.cc

* HAVE_JN_ ?

* TODO: Update the p4est-setup.sh script...

