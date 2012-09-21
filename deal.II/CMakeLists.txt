SET(VERSION "8.0.pre")

#
# This is the main CMakeLists.txt file for the deal.II project
#
# The main purpose of this file is to set all the configuration options and
# then call a lot of modules residing under ./cmake/*/ that do the
# hard work.
#
#
# Note: All further explained options can be either set via ccmake or at
# the promp via:
#
#   $ cmake -D<OPTION>=<VALUE> [...] path/to/source/tree
#
#
# The first bunch of options of type DEAL_II_WITH_<feature> are for
# configuration of the various features of the deal.II library.
# An enabled feature will usually add one ore more dependencies to external
# or 'contrib' (i.e. bundled with deal.II and residing under ./contrib/)
# libraries.
#
#
# The option DEAL_II_ALLOW_CONTRIB determines the behaviour of the
# dependency resolution:
#
# - If set to ON external libraries still have precedence. But if there is
#   no external library the bundled contrib library will be used.
#
# - If set to OFF the contrib libraries will not be used and the dependency
#   resolution will fail if there is no external library. EXCEPT (there is
#   always an exception, isn't it?) DEAL_II_FORCE_CONTRIB_<library> forces
#   the use of the bundled contrib library.
#
#
# Autoconfiguration: If the option DEAL_II_FEATURE_AUTODETECTION is enabled
# the DEAL_II_WITH_<feature> toggles will be automatically set (overwriting
# any previous configuration) depending on whether they can be supported or
# not.
# (Note: DEAL_II_FEATURE_AUTODETECTION will respect DEAL_II_ALLOW_CONTRIB
# and DEAL_II_FORCE_CONTRIB_<library>)
#
#
# The option DEAL_II_FORCE_CONTRIB_<library> forces the use of the bundled
# contrib library regardless whether DEAL_II_ALLOW_CONTRIB is set to OFF or
# an external library is found.
#
# REMARKS:
#
# - To ensure that no bundled library is used at all DEAL_II_ALLOW_CONTRIB,
#   as well as every DEAL_II_FORCE_CONTRIB_<library> have to be set to OFF.
#
# - DEAL_II_FORCE_CONTRIB_<library>=ON will _not_ automatically enable the
#   corresponding DEAL_II_WITH_<feature> toggle. This has to be set
#   separately.
#
#
# And, last but not least, CMAKE_BUILD_TYPES:
#
# We support the "Debug" and "Release" build targets. Default is the
# "Release" target.
#
# - Debug will enable a lot of debug code paths and Assertions, as well as
#   compiling the library with a suitable set of debug compiler flags.
#
# - Release will chose a suitable set of optimizing flags.
#
#
# deal.II will configure sensible default CFLAGS and CXXFLAGS depending on
# platform, compiler and build target. There are two options, if this is
# not desired:
#
# - Option 1: Disable the configuration completely by setting
#   DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS to OFF.
#   Beware of the fact that certain features may still pull in necessary
#   compiler flags.
#
# - Option 2: Overwrite the default configuration by setting the following
#   cached variables:
#
#     CMAKE_CXX_FLAGS         -  used during all builds
#     CMAKE_CXX_FLAGS_DEBUG   -  additional flags used during debug builds
#     CMAKE_CXX_FLAGS_RELEASE -  additional flags used during release builds
#     (same for CMAKE_C_...)
#
# The content of the cached variables will be preserved and added
# *_TO THE END_* of the default compiler flags, hence giving a possibility
# for overwriting a flag. E.g.: -Wsign-compare, set by the build system,
# can be overwritten by specifying:
#
#   $ cmake -DCMAKE_CXX_FLAGS="-Wno-sign-compare" <...>
#
#
# Build options:
#
# - BUILD_SHARED_LIBS and
# - DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS
#
# Installation:
#
# - DEAL_II_INSTALL_EXAMPLES
# - DEAL_II_INSTALL_COMPAT_FILES
# - DEAL_II_INSTALL_PROJECT_CONFIG
#
# - Explain Components and Options...
# - Explain how to build and install only a specific component:
#
#   $ make <component>
#   $ cmake -DCOMPONENT=<component> -P build/cmake_install.cmake
#




###########################################################################
#                                                                         #
#                     General configuration options:                      #
#                                                                         #
###########################################################################

#
# Feature selection: external libraries and where to take them (if there
# is a choice, e.g. from contrib vs external)
#

# TODO: Some of the description strings span over multiple lines. Fix
# this somehow.

OPTION(DEAL_II_FEATURE_AUTODETECTION
  "Enables feature autodetection. This will automatically overwrite all
  DEAL_II_WITH_<...> toggles depending on whether they can be supported or
  not."
  ON)

OPTION(DEAL_II_ALLOW_CONTRIB
  "Allow the use of contrib libraries bundled with the source tarball.
  Note: If set to off DEAL_II_FORCE_CONTRIB* will still pull in bundled
  packages, so to ensure that only external libraries are used
  DEAL_II_ALLOW_CONTRIB as well as all DEAL_II_FORCE_CONTRIB_* have to be
  OFF"
  OFF)

OPTION(DEAL_II_WITH_ARPACK
  "Build deal.II with support for ARPACK."
  OFF)

OPTION(DEAL_II_WITH_BLAS
  "Build deal.II with support for BLAS."
  OFF)

SET(DEAL_II_WITH_BOOST "ON"
  CACHE BOOL "Build deal.II with support for boost." FORCE
  ) # Always true. We need it :-]
OPTION(DEAL_II_FORCE_CONTRIB_BOOST
  "Always use the bundled boost library instead of an external one."
  ON)

OPTION(DEAL_II_WITH_DOXYGEN
  "Enable the generation and installation of the documentation via doxygen.
  This feature adds a COMPONENT \"documentation\" to the build system."
  OFF)

OPTION(DEAL_II_WITH_FUNCTIONPARSER
  "Build deal.II with support for functionparser."
  OFF)
OPTION(DEAL_II_FORCE_CONTRIB_FUNCTIONPARSER
  "Always use the bundled functionparser library instead of an external one."
  ON)

OPTION(DEAL_II_WITH_LAPACK
  "Build deal.II with support for LAPACK."
  OFF)

OPTION(DEAL_II_WITH_METIS
  "Build deal.II with support for Metis."
  OFF)

OPTION(DEAL_II_WITH_MPI
  "Build deal.II with support for mpi."
  OFF)

OPTION(DEAL_II_WITH_NETCDF
  "Build deal.II with support for netcdf."
  OFF)

OPTION(DEAL_II_WITH_P4EST
  "Build deal.II with support for P4EST."
  OFF)

OPTION(DEAL_II_WITH_TBB
  "Build deal.II with support for tbb. This will enable thread support in deal.II."
  OFF)
OPTION(DEAL_II_FORCE_CONTRIB_TBB
  "Always use the bundled tbb library instead of an external one."
  OFF)

OPTION(DEAL_II_WITH_TRILINOS
  "Build deal.II with support for trilinos."
  OFF)

OPTION(DEAL_II_WITH_UMFPACK
  "Build deal.II with support for UMFPACK."
  OFF)
OPTION(DEAL_II_FORCE_CONTRIB_UMFPACK
  "Always use the bundled umfpack library instead of an external one."
  OFF)

OPTION(DEAL_II_WITH_ZLIB
  "Build deal.II with support for zlib."
  OFF)

#
# Build configuration: configuration options regarding compilation and
# installation of the deal.II library
#

OPTION(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS
  "configure sensible default CFLAGS and CXXFLAGS depending on platform,
  compiler and build target. Disabling this option will turn off almost all
  compiler flag setup."
  ON)

OPTION(DEAL_II_INSTALL_EXAMPLES
  "Enable the installation of the example steps.
  This feature adds a COMPONENT \"examples\" to the build system."
  ON)

OPTION(DEAL_II_INSTALL_COMPAT_FILES
  "Enable the installation of the example steps.
  This feature adds a COMPONENT \"compat_files\" to the build system."
  ON)

OPTION(DEAL_II_INSTALL_PROJECT_CONFIG
  "Enable the installation of a cmake project config. Useful for
  autodetection ( FIND_PACKAGE(DEAL_II) ) of the deal.II library out of
  cmake build systems.
  This feature adds a COMPONENT \"project_config\" to the build system."
  ON)




###########################################################################
#                                                                         #
#                             Configuration:                              #
#                                                                         #
###########################################################################


#
# General configuration for cmake:
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

SET(CMAKE_INCLUDE_CURRENT_DIR true)

SET(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  ${CMAKE_SOURCE_DIR}/cmake/
  ${CMAKE_SOURCE_DIR}/cmake/modules/
  )

#
# Load all macros:
#
FILE(GLOB macro_files "cmake/macros/*.cmake")
FOREACH(file ${macro_files})
  MESSAGE(STATUS "Include ${file}")
  INCLUDE(${file})
ENDFOREACH()

INCLUDE(setup_external_macros)


#
# We have to initialize some cached variables before PROJECT is called, so
# do it at this point:
#
INCLUDE(setup_cached_variables)


#
# Now, set the project and setup the rest:
#
PROJECT(deal.II)

INCLUDE(setup_deal_ii)

INCLUDE(setup_compiler_flags)


##########################################################################
#
#                             FAT NOTE:
#
# To keep things clean, only the following cmake variables should be
# altered in the platform checks and feature configuration
# (beside setting a lot of DEAL_II_* definitions...):
#
#
# For internal use:
#
#     CMAKE_SHARED_LINKER_FLAGS
#         for setting necessary linker flags for the deal.II library.
#
#     CMAKE_CXX_FLAGS (and CMAKE_C_FLAGS)
#         for setting necessary compiler flags, e.g. -std=c++11 (if
#         available).
#
#     INCLUDE_DIRECTORIES()
#         For setting necessary include dirs for the compilation of the
#         deal.II library.
#
#     ADD_DEFINITIONS()
#         For setting necessary preprocessor definitions ("-D<...>") for
#         the compilation of the deal.II library
#
ADD_CUSTOM_TARGET(deal_ii_target_dependencies)
#         Used as a dummy target for all the stuff that has to be done
#         before the library can be compiled.
#
#
# For internal and external use:
#
#     DEAL_II_EXTERNAL_LIBRARIES
#         Used to keep track of external libraries, the deal.II library
#         and user programs have to be linked against.
#
#
# For external use:
#
#     DEAL_II_EXTERNAL_DEFINITIONS
#         Used to keep track of external preprocessor definitions,
#         necessary for the compilation of user programs.
#
#
#     DEAL_II_EXTERNAL_INCLUDE_DIRS
#         Used to keep track of external include dirs, necessary for the
#         compilation of user programs.
#
##########################################################################


#
# Run all system checks:
#
FILE(GLOB check_files "cmake/check/*.cmake")
FOREACH(file ${check_files})
  MESSAGE(STATUS "Include ${file}")
  INCLUDE(${file})
ENDFOREACH()


#
# Feature configuration:
#
FILE(GLOB configure_files "cmake/configure/configure_*.cmake")
# Make sure the configure_*.cmake files are included in alphabetical order
# to fullfill some dependencies...
LIST(SORT configure_files)
FOREACH(file ${configure_files})
  MESSAGE(STATUS "Include ${file}")
  INCLUDE(${file})
ENDFOREACH()


#
# Finalize the configuration:
#
INCLUDE(setup_finalize)




###########################################################################
#                                                                         #
#                    Compile and install the library:                     #
#                                                                         #
###########################################################################


#
# Build the toolchain and prepare the includes:
#
ADD_SUBDIRECTORY(config)
ADD_SUBDIRECTORY(scripts)
ADD_SUBDIRECTORY(include)

#
# Compile and install the library:
#
ADD_SUBDIRECTORY(source)


#
# Prepare some directories for installation
#
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(examples)

#
# Define some custom targets for convenience
#
INCLUDE(setup_custom_targets)
