Some last bits and thoughts:

* What about migrating the testsuite to ctest/cdash? This would have
  several advantages
    - good integration with the cmake build system (*cough*)
    - automatic build (!) and regression tests
    - cdash is quite nice and we don't have to write the web front end by
      hand...

* A TECPLOT feature has to be written

* Enable compilation of bundled TBB under Windows/Cygwin and Windows/MinGW

* Clean up Platform introspection (search for sse2 and avx) and document
  how to disable it (important for building for a different target system)

* Provide a configuration template for cross compiling with CMake,

* Set VERSION in CMakeLists.txt or in an extra file?

* Some failing tests:

    2012-12-05 19:00  2   base/thread_validity_08
    2012-12-05 19:00  2   base/thread_validity_09
    2012-12-05 19:01  3   deal.II/dof_renumbering_04
    2012-12-05 19:01  3   deal.II/dof_renumbering_04a
    2012-12-05 19:01  3   deal.II/support_point_map

        * This stuff looks bad...

    2012-12-05 19:00  3   base/utilities_01

        * Random number generator

    2012-12-05 19:00  0   base/path_search

        * broken because DEAL_II_PATH is not available any more in config.h


