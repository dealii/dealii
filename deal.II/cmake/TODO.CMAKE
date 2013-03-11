Major:

* Mumps may have a transitive link interface to Scotch. Fix this

* Get rid of the last bits of unguarded platform introspection

Some last bits and thoughts:

* A TECPLOT feature has to be written

* Enable compilation of bundled TBB under Windows/Cygwin and Windows/MinGW

* Provide a configuration template for cross compiling with CMake,

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

* document that we no longer need the two-lib installation of p4est

* adjust doc/development/testsuite.html

* update https://code.google.com/p/dealii/wiki/Windows

