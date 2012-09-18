SET(VERSION "8.0.pre")

#
# This is the main CMakeLists.txt file for the deal.II project
#
# The main purpose of this file is to set all the configuration toggles and
# then call a lot of modules residing under ./contrib/cmake/*/ that do the
# hard work.
#
# The first bunch of options of type DEAL_II_WITH_<feature> are for
# configuration of the various features of the deal.II library.
# An enabled feature will usually add one ore more dependencies to external
# or 'contrib' (i.e. bundled with deal.II and residing under ./contrib/)
# libraries.
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
# We support the "Debug" and "Release" build types. Default is the
# "Release" target.
#
# Debug will enable a lot of debug code paths and Assertions, as well as
# compiling the library with a suitable set of debug compiler flags.
#
# Release will chose a suitable set of optimizing flags.
#
# The default CXXFLAGS (and CFLAGS) are stored in the following cached
# variables (as is common practise for cmake):
#
#   CMAKE_CXX_FLAGS         -  used during all builds
#   CMAKE_CXX_FLAGS_DEBUG   -  additional flags used during debug builds
#   CMAKE_CXX_FLAGS_RELEASE -  additional flags used during release builds
#
#   (same for CMAKE_C_...)
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

OPTION(DEAL_II_FEATURE_AUTODETECTION
  "Enables feature autodetection. This will automatically overwrite all
  DEAL_II_WITH_<...> toggles depending on whether they can be supported or
  not."
  OFF)

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

OPTION(DEAL_II_WITH_FUNCTIONPARSER
  "Build deal.II with support for functionparser."
  OFF)
OPTION(DEAL_II_FORCE_CONTRIB_FUNCTIONPARSER
  "Always use the bundled functionparser library instead of an external one."
  OFF)

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

OPTION(DEAL_II_FORCE_CONTRIB_BOOST
  "Always use the bundled boost library instead of an external one."
  OFF)




###########################################################################
#                                                                         #
#                             Configuration:                              #
#                                                                         #
###########################################################################


#
# General configuration for cmake:
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

SET(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/contrib/cmake/"
  "${CMAKE_SOURCE_DIR}/contrib/cmake/modules/"
  )

#
# Load all macros:
#
FILE(GLOB macro_files "contrib/cmake/macros/*.cmake")
FOREACH(file ${macro_files})
  MESSAGE(STATUS "Include ${file}")
  INCLUDE(${file})
ENDFOREACH()

INCLUDE(setup_external_macros)


#
# We have to initialize some cached variable before PROJECT is called, so
# configure the cached, user editable CMAKE_BUILD_TYPES, CMAKE_CXX_FLAGS,
# etc. variables at this point:
#
INCLUDE(setup_cached_variables)

#
# Now, set the project and configure some deal_II specific variables:
#
PROJECT(deal.II)
INCLUDE(setup_deal_ii)


#
#   ####################
#   #     FAT NOTE:    #
#   ####################
#
# To keep things clean, only the following cmake variables should be
# altered in the platform checks and features configuration
# (beside setting a lot of DEAL_II_* definitions...):
#
# CMAKE_SHARED_LINKER_FLAGS
#
#   for setting necessary linker flags for the deal.II library.
#
# CMAKE_C_FLAGS and CMAKE_CXX_FLAGS (uncached)
#
#   for setting necessary compiler flags, e.g. -std=c++11 (if available).
#   Note: Uncached variables are not visible to the user, and cannot be
#   altered.
#
ADD_CUSTOM_TARGET(deal_ii_target_dependencies)
#
#   Used as a dummy target for all the stuff that has to be done before the
#   library can be compiled.
#
SET(deal_ii_external_libraries)
#
#   Used to keep track of external libraries, against which the deal.II
#   library will be linked.
#


#
# Run all system checks:
#

FILE(GLOB check_files "contrib/cmake/check/*.cmake")
FOREACH(file ${check_files})
  MESSAGE(STATUS "Include ${file}")
  INCLUDE(${file})
ENDFOREACH()


#
# Feature configuration:
#
FILE(GLOB configure_files "contrib/cmake/configure/configure_*.cmake")

# Make sure the configure_*.cmake files are included in alphabetical order:
LIST(SORT configure_files)

FOREACH(file ${configure_files})
  MESSAGE(STATUS "Include ${file}")
  INCLUDE(${file})
ENDFOREACH()




###########################################################################
#                                                                         #
#                          Compile the library:                           #
#                                                                         #
###########################################################################


#
# Build the toolchain:
#

ADD_SUBDIRECTORY(common)

# TODO
# CONFIGURE_FILE(
#   ${CMAKE_CURRENT_SOURCE_DIR}/Version.in
#   ${CMAKE_CURRENT_BINARY_DIR}/Version
#   )


#
# Compile the library:
#

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_BINARY_DIR}/include/
  )

ADD_SUBDIRECTORY(include)

ADD_SUBDIRECTORY(source)

MESSAGE(STATUS "Configured features:")
GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  IF(var MATCHES "DEAL_II_WITH")
    MESSAGE(STATUS "${var} : ${${var}}")
  ENDIF()
ENDFOREACH()
