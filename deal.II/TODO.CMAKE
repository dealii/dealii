Major:

* Furnish the MPI dependencies (arpack, hdf5, petsc, p4est, trilinos, metis)

* Add the rest of the Bug tests

* Update the online documentation on

    - How to install deal.II

    - How to interface with external libraries (Unify this...)

    - How to use deal.II
      * Autopilot CMakeLists.txt
      * minimal CMakeLists.txt
      * direct interface

    - The documentation in the example steps

* Port the testsuite to cmake


Minor:

* The following features have to be written:
    * TECPLOT

* (Maybe) re-enable validate-xrefs

* Fixup the TODOs in
    source/lac/sparse_dircect.cc
    source/base/path_find.cc
    source/base/function_parser.cc

* HAVE_JN_ ?

* TODO: Update the p4est-setup.sh script...

