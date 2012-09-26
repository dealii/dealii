


          ==================================================

            The cmake build system for the deal.II project

          ==================================================

          Author: Matthias Maier, University Heidelberg 2012




=================
TABLE OF CONTENTS
=================

  Introduction

      The philosophy

      Quick start

      Useful bits

  Configuration

      Feature configuration

      External library locations

      Component selection

      Build configuration

      Installation

  How to use deal.II in your cmake project

      Finding the deal.II library

      Setting up necessary configuration variables

  Development

      Fat note regarding platform checks and feature configuration

      Writing platform checks

      Writing features




============
INTRODUCTION
============


  The philosophy
  ==============

          »I need a new keyboard, my CAPS-KEY is stuck because of writing all
           this cmake script...«


  Quick start
  ===========


  * Make a build directory, configure, compile and install:

        $ mkdir build
        $ cd build
        $ cmake -DCMAKE_INSTALL_PREFIX=/path/install/dir ../
        $ make install

  * Check out the examples, shipped with deal.II:

        $ cd /path/to/install/dir/examples


  Useful bits
  ===========


  * All cached variables can be either set via

      ( $ cd build )
        $ ccmake

    or at the prompt via:

      ( $ cd build )
        $ cmake -D<OPTION>=<VALUE> [...] path/to/source/tree

  * Configure with autodetection, afterwards disable a feature:

      ( $ cd build )
        $ cmake -DDEAL_II_FEATURE_AUTODETECTION=ON ../
        $ cmake -DDEAL_II_FEATURE_AUTODETECTION=OFF -DDEAL_II_WITH_METIS=OFF ../
        $ make install

  * Build and install only a specific component:

      ( $ cd build )
        $ make <component>
        $ cmake -DCOMPONENT=<component> -P build/cmake_install.cmake

  * The build system cmake generates allows fine grained control. It is
    e.g. possible to only compile a single file:

        $ cd source/dofs
        $ make obj_dofs.inst
        [  0%] Building CXX object scripts/CMakeFiles/expand_instantiations.dir/expand_instantiations.cc.o
        Linking CXX executable expand_instantiations
        [  0%] Built target expand_instantiations
        [ 33%] Generating block_info.inst
        [ 33%] Generating dof_accessor.inst
        [ 33%] Generating dof_handler.inst
        [ 66%] Generating dof_handler_policy.inst
        [ 66%] Generating dof_objects.inst
        [ 66%] Generating dof_renumbering.inst
        [100%] Generating dof_tools.inst
        [100%] Built target obj_dofs.inst
        $ make dof_handler.o
        Building CXX object source/dofs/CMakeFiles/obj_dofs.dir/dof_handler.cc.o

    For a list of valid targets of a given folder just type

        $ make help




=============
CONFIGURATION
=============


  The various configuration options of the deal.II library are organized in
  three categories: Feature, Component and build/install configuration.


  Feature configuration
  =====================

  deal.II provides (optional) interfaces to quite a number of external
  libraries:

    DEAL_II_WITH_ARPACK
        Build deal.II with support for ARPACK

    DEAL_II_WITH_BLAS
        Build deal.II with support for BLAS

    DEAL_II_WITH_BOOST
        Build deal.II with support for boost. Always activated

    DEAL_II_WITH_DOXYGEN
        Build deal.II with support for doxygen and dot

    DEAL_II_WITH_FUNCTIONPARSER
        Build deal.II with support for functionparser

    DEAL_II_WITH_LAPACK
        Build deal.II with support for LAPACK

    DEAL_II_WITH_METIS
        Build deal.II with support for Metis

    DEAL_II_WITH_MPI
        Build deal.II with support for mpi

    DEAL_II_WITH_NETCDF
        Build deal.II with support for netcdf

    DEAL_II_WITH_P4EST
        Build deal.II with support for p4est

    DEAL_II_WITH_PETSC
        Build deal.II with support for petsc

    DEAL_II_WITH_TBB
        Build deal.II with support for tbb. This will enable thread support in deal.II

    DEAL_II_WITH_TRILINOS
        Build deal.II with support for trilinos

    DEAL_II_WITH_UMFPACK
        Build deal.II with support for UMFPACK

    DEAL_II_WITH_ZLIB
        Build deal.II with support for zlib

  An enabled feature will usually add one ore more dependencies to external
  or 'contrib' (i.e. bundled with deal.II and residing under ./contrib/)
  libraries.


  There are some options to determine the behaviour of the dependency
  resolution:

  The option DEAL_II_ALLOW_CONTRIB:

    - If set to ON external libraries still have precedence. But if there is
      no external library the bundled contrib library will be used.

    - If set to OFF the contrib libraries will not be used and the dependency
      resolution will fail if there is no external library. EXCEPT (there is
      always an exception, isn't it?) DEAL_II_FORCE_CONTRIB_<library> forces
      the use of the bundled contrib library.

  The option DEAL_II_FORCE_CONTRIB_(BOOST|FUNCTIONPARSER|TBB|UMFPACK)
  forces the use of the bundled  contrib library regardless whether
  DEAL_II_ALLOW_CONTRIB is set to OFF or an external library is found.


  Autoconfiguration:

  If the option DEAL_II_FEATURE_AUTODETECTION is enabled the
  DEAL_II_WITH_<feature> toggles will be automatically set (overwriting any
  previous configuration) depending on whether they can be supported or
  not.
  (Note: DEAL_II_FEATURE_AUTODETECTION will respect DEAL_II_ALLOW_CONTRIB
  and DEAL_II_FORCE_CONTRIB_<feature>)

  REMARKS:

  - To ensure that no bundled library is used at all DEAL_II_ALLOW_CONTRIB,
    as well as every DEAL_II_FORCE_CONTRIB_<library> have to be set to OFF.

  - DEAL_II_FORCE_CONTRIB_<library>=ON will _not_ automatically enable the
    corresponding DEAL_II_WITH_<feature> toggle. This has to be set
    separately.


  External library locations
  ==========================

  External libraries will be searched depending on hints in the following
  order:

  1.) Paths specified via CMAKE_PREFIX_PATH take precedence, e.g. with

          $ cmake -DCMAKE_PREFIX_PATH=~/workspace/local ../deal.II

      libraries from ~/workspace/local will be preferred for dependency
      resolution.

  2.) Hints given by <library>_DIR via command line or environment for
      _some_ libraries:

          $ cmake -DP4EST_DIR=~/workspace/p4est-install/ ../deal.II
      or
          $ export P4EST_DIR=~/workspace/p4est-install/
          $ cmake ../deal.II
      where -D<library>_DIR takes precedence over environment.

      Currently, the following variables will be considered:

        METIS_DIR,
        P4EST_DIR (and SC_DIR),
        PETSC_DIR and PETSC_ARCH (forming ${PETSC_DIR}/${PETSC_ARCH}),
        TRILINOS_DIR,
        UMFPACK_DIR (and AMD_DIR)

  3.) The system default locations for libraries and includes.


  Alternatively, cached variables set by the Find<Module> mechanism may be
  set, hinted or overwritten directly!

  These variables are usually called <library>_INCLUDE_DIR(|S) and
  <library>_(LIBRARY|LIBRARIES) (highly dependend of the actual library...).

  You can get a list via

      $ ccmake ../deal.II (or $ make edit_cache)

  and entering advanced configuration mode by pressing [t].
  Variables that could not be determined are suffixed with -NOTFOUND and
  may be set by hand.


  Component selection
  ===================


  The following options control which components of deal.II will be
  configured, build and installed:

    DEAL_II_COMPONENT_EXAMPLES:
      Enable configuration and installation of the example steps.
      This adds a COMPONENT "examples" to the build system.

    DEAL_II_COMPONENT_DOCUMENTATION:
      Enable configuration, build and installation of the documentation.
      This adds a COMPONENT "documentation" to the build system.

    DEAL_II_COMPONENT_COMPAT_FILES:
      Enable installation of legacy files and tools for compatibility with
      the old build system
      This adds a COMPONENT "compat_files" to the build system.

    DEAL_II_COMPONENT_PROJECT_CONFIG:
      Enable configuration and installation of a cmake project config.
      Useful for finding ( FIND_PACKAGE(DEAL_II) ) the deal.II library out
      of cmake build systems. This adds a COMPONENT "project_config" to the
      build system.




  Build configuration
  ===================


  CMAKE_BUILD_TYPE: We support the "Debug", "Release" and "DebugRelease"
  build targets. Default is the "DebugRelease" target.

  - Debug will compile a debug library libdeal_II.g.so with a lot of debug
    code paths and assertions enabled, as well as a suitable set of debug
    compiler flags.

  - Release will build an optimized libdeal_II.so.

  - DebugRelease will build both libraries.


  deal.II will configure sensible default CFLAGS and CXXFLAGS depending on
  platform, compiler and build target. There are two options, if this is
  not desired:

  - Option 1: Disable the configuration completely by setting
    DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS to OFF.
    Beware of the fact that certain features may still pull in necessary
    compiler flags.

  - Option 2: Overwrite the default configuration by setting the following
    cached variables:

      CMAKE_CXX_FLAGS           - used during all builds
      DEAL_II_CXX_FLAGS_DEBUG   - additional flags for the debug library
      DEAL_II_CXX_FLAGS_RELEASE - additional flags for the release library

  The content of the cached variables will be preserved and added
  *_TO THE END_* of the default compiler flags, hence giving a possibility
  for overwriting a flag. E.g.: -Wsign-compare, set by the build system,
  can be overwritten by specifying:

    $ cmake -DCMAKE_CXX_FLAGS="-Wno-sign-compare" <...>


  The build can be further controlled by the following variables:

  - BUILD_SHARED_LIBS: If set, deal.II will be linked as a shared library

  - CMAKE_INSTALL_RPATH_USE_LINK_PATH: If set, the deal.II library will be
    installed with rpaths  set for all libraries outside of the system
    search paths




  Installation
  ============


  CMAKE_INSTALL_PREFIX determines the location, where the deal.II library
  will be installed.

  Please note that depending on whether DEAL_II_COMPONENT_COMPAT_FILES is
  set, there will be different folder structures:

  With DEAL_II_COMPONENT_COMPAT_FILES=ON:

  ${CMAKE_INSTALL_PREFIX}/
      common
      common/scripts
      doc
      examples
      include
      lib


  With DEAL_II_COMPONENT_COMPAT_FILES=OFF:

  ${CMAKE_INSTALL_PREFIX}/
      include
      lib${LIB_SUFFIX}
      /share/doc/deal.II/examples
      /share/doc/deal.II/html


  The individual target directories can be overwritten by
      DEAL_II_CMAKE_MACROS_RELDIR
      DEAL_II_DOCUMENTATION_RELDIR
      DEAL_II_EXAMPLES_RELDIR
      DEAL_II_INCLUDE_RELDIR
      DEAL_II_LIBRARY_RELDIR
      DEAL_II_PROJECT_CONFIG_RELDIR




========================================
How to use deal.II in your cmake project
========================================


  Finding the deal.II library
  ===========================


  If DEAL_II_COMPONENT_PROJECT_CONFIG was set, finding the deal.II library
  should be no more than

    FIND_PACKAGE(deal.II CONFIG REQUIRED)

  in your CMakeLists.txt file.

  You may have to hint for the location of the deal.IIConfig.cmake file.
  Either by directly specifying deal.II_DIR to point to the path were the
  deal.IIConfig.cmake files are located:

    $ cmake -Ddeal.II_DIR=/path/to/the/config/file <...>

  Or by specifying a search path via CMAKE_PREFIX_PATH, e.g.

    $ cmake -DCMAKE_PREFIX_PATH=~/workspace/local

  In this case deal.IIConfig.cmake will be searched in
    ~/workspace/local/
    ~/workspace/local/lib/cmake/deal.II/


  Setting up necessary configuration variables
  ============================================

  Importing the deal.IICONFIG.cmake file via FIND_PACKAGE will set a bunch
  of variables. All of the form DEAL_II_*.

  Usually for actually using deal.II the following configuration steps are
  necessary. This can be either done by hand (a), or set up via macros (b).


  a) Configuration by hand:

    - Set up all include directories for headers:

        INCLUDE_DIRECTORIES(
          # The deal.II include directories:
          ${DEAL_II_INCLUDE_DIRS}
          # External include directories necessary for user programs:
          ${DEAL_II_EXTERNAL_INCLUDE_DIRS}
          )

    - deal.II usually ships with an optimized Release and a Debug version
      of the library. So it is a good idea to set up CMAKE_BUILD_TYPE
      accordingly:

        SET(CMAKE_BUILD_TYPE "Debug" CACHE
          "Choose the type of build, options are: Debug, Release"
          )

      A cached variable ensures that we can later switch the build type
      by editing the cache:

        $ ccmake <...> or $ make edit_cache

    - Often, it is a good idea to use the same compiler and linker
      flags as the deal.II library.

        SET(CMAKE_CXX_FLAGS ${DEAL_II_CXX_FLAGS})
        SET(CMAKE_CXX_FLAGS_RELEASE ${DEAL_II_CXX_FLAGS_RELEASE})
        SET(CMAKE_CXX_FLAGS_DEBUG ${DEAL_II_CXX_FLAGS_DEBUG})

      (Optionally you can set them up with 'CACHE "..."' to be able
      to edit them via ccmake or make edit_cache)

    - _After_ this, specify your project name:

        PROJECT(myProject)

      This ensures that the compiler detection and platform setup that
      is issued by calling PROJECT(...) will run after we have set up
      our cached variables. This way it is our choice of variables that
      will be set and not the one determined in PROJECT(...)will be set
      and not the one determined in PROJECT(...)

    - After defining your targets, e.g.

        ADD_EXECUTABLE(step-1 step-1.cc)

      we have to specify to link against deal.II and all external
      libraries:

        TARGET_LINK_LIBRARIES(step-1
          debug
            ${DEAL_II_LIBRARIES_DEBUG}
            ${DEAL_II_EXTERNAL_LIBRARIES_DEBUG}
          optimized
            ${DEAL_II_LIBRARIES_RELEASE}
            ${DEAL_II_EXTERNAL_LIBRARIES_RELEASE}
          general
            ${DEAL_II_EXTERNAL_LIBRARIES}
          )

     - As well as using some preprocessor definitions:

        SET_TARGET_PROPERTIES(step-1 PROPERTIES
          COMPILE_DEFINITIONS
            "${DEAL_II_USER_DEFINITIONS}"
          COMPILE_DEFINITIONS_DEBUG
            "${DEAL_II_USER_DEFINITIONS_DEBUG}"
          COMPILE_DEFINITIONS_RELEASE
            "${DEAL_II_USER_DEFINITIONS_RELEASE}"
          )

     - Optionally, you can set the link flags to the one used by the
       deal.II library:

         SET_TARGET_PROPERTIES(step-1 PROPERTIES
           LINK_FLAGS
             "${DEAL_II_LINKER_FLAGS}"
           LINK_FLAGS_DEBUG
             "${DEAL_II_LINKER_FLAGS_DEBUG}"
           LINK_FLAGS_RELEASE
             "${DEAL_II_LINKER_FLAGS_RELEASE}"
           )

      And this it is.


  b) Configuration with the help of two convencience macros:

    All the steps explained in (a) can be automatically set with the help
    of two convenience macros. This boils down to the following example
    code:

        FIND_PACKAGE(deal.II CONFIG REQUIRED)

        DEAL_II_INITIALIZE_CACHED_VARIABLES()

        PROJECT(step-1)

        ADD_EXECUTABLE(step-1 step-1.cc)
        DEAL_II_SETUP_TARGET(step-1)




===========
DEVELOPMENT
===========


  Fat note regarding platform checks and feature configuration
  ============================================================


                            ============

                              FAT NOTE

                            ============


  To keep things clean, only the following variables should appended in the
  platform checks and feature configuration (beside setting a lot of
  DEAL_II_* definitions...):


  For internal use:

      CMAKE_SHARED_LINKER_FLAGS
      DEAL_II_SHARED_LINKER_FLAGS_DEBUG
      DEAL_II_SHARED_LINKER_FLAGS_RELEASE

          for setting necessary linker flags for the deal.II library.

      CMAKE_CXX_FLAGS (and CMAKE_C_FLAGS)
      DEAL_II_CXX_FLAGS_DEBUG (and DEAL_II_C_FLAGS_DEBUG)
      DEAL_II_CXX_FLAGS_RELEASE (and DEAL_II_C_FLAGS_RELEASE)

          for setting necessary compiler flags, e.g. -std=c++11 (if
          available).

      INCLUDE_DIRECTORIES(...)

          For setting necessary include dirs for the compilation of the
          deal.II library.

      DEAL_II_DEFINITIONS
      DEAL_II_DEFINITIONS_DEBUG
      DEAL_II_DEFINITIONS_RELEASE

          For setting necessary preprocessor definitions ("-D<...>") for
          the compilation of the deal.II library


  For internal and external use:

      DEAL_II_EXTERNAL_LIBRARIES
      DEAL_II_EXTERNAL_LIBRARIES_DEBUG
      DEAL_II_EXTERNAL_LIBRARIES_RELEASE

          Used to keep track of external libraries, the deal.II library
          and user programs have to be linked against.


  For external use:

      DEAL_II_USER_DEFINITIONS
      DEAL_II_USER_DEFINITIONS_DEBUG
      DEAL_II_USER_DEFINITIONS_RELEASE

          Used to keep track of external preprocessor definitions,
          necessary for the compilation of user programs.


      DEAL_II_USER_INCLUDE_DIRS

          Used to keep track of external include dirs, necessary for the
          compilation of user programs.




  Writing platform checks
  =======================


  Platform checks reside under cmake/check.
  We currently have

    cmake/check/check_for_compiler_bugs.cmake

    cmake/check/check_for_compiler_features.cmake

    cmake/check/check_for_cxx_features.cmake


  General notes:

  - A platform check should have a prominent comment explaining what
    it does and why it is there, as well as a stating the author and the
    year.

  - Definition toggles in include/deal.II/base/config.h.in should have
    a prominent comment explaining it and should be grouped by file
    exporting the definition


  Writing features
  ================

  For a feature <feature> the following the following toggles, variables and
  macros should be defined:

  (Note: <FEATURE> means all caps, <feature> means all lowercase)


  In the top level ./CMakeLists.txt:

    DEAL_II_WITH_<FEATURE> (bool, mandatory)
      An option for enabling or disabling the configuration of <feature>

    DEAL_II_FORCE_CONTRIB_<FEATURE> (bool, optional)
      If <feature> can be set up by contrib dependencies, this
      configuration option must be present to force the use of the contrib
      dependencies


  A file cmake/configure/configure_<feature>.cmake defining how to
  configure <feature>:

    DEAL_II_WITH_<FEATURE> (bool, mandatory)
       If DEAL_II_FEATURE_AUTODETECTION is OFF, this boolean determines
       whether the feature will be configured.
       If DEAL_II_FEATURE_AUTODETECTION is ON, this boolean will
       automatically be set if configuring the feature was successful.

    FEATURE_<FEATURE>_DEPENDS (variable, optional)
       a variable which contains an optional list of other features
       this feature depends on (and which have to be enbled for this feature
       to work.) The features must be given with the full option toggle:
       DEAL_II_WITH_[...]

    FEATURE_<FEATURE>_HAVE_CONTRIB  (variable, optional)
       which should either be set to TRUE if all necessary libraries of the
       features comes bundled with deal.II and hence can be supported
       without external dependencies, or unset.

    FEATURE_<FEATURE>_CONFIGURE_CONTRIB(var)  (macro, optional)
       which should setup all necessary configuration for the feature with
       contrib source dependencies. var set to TRUE indicates success,
       otherwise this script gives an error.

    FEATURE_<FEATURE>_FIND_EXTERNAL(var)  (macro, mandatory)
       which should set var to TRUE if all dependencies for the feature are
       fulfilled. In this case all necessary variables for
       FEATURE_<FEATURE>_CONFIGURE_EXTERNAL must be set. Otherwise
       var should remain unset.
       This macro must not give an error (SEND_ERROR or FATAL_ERROR).

    FEATURE_<FEATURE>_CONFIGURE_EXTERNAL(var)  (macro, mandatory)
       which should setup all necessary configuration for the feature with
       external dependencies. var set to TRUE indicates success,
       otherwise this script gives an error.

    FEATURE_<FEATURE>_CUSTOM_ERROR_MESSAGE() (variable, optional)
       which should either be set to TRUE if FEATURE_<FEATURE>_ERROR_MESSAGE
       is set up, or be undefined.

    FEATURE_<FEATURE>_ERROR_MESSAGE()  (macro, optional)
       which should print a meaningful error message (with SEND_ERROR) for
       the case that no external library was found (and contrib is not
       allowed to be used.) If not defined, a suitable default error message
       will be printed.


