PROJECT(deal.II)

SET(VERSION "8.0.pre")

#
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

#
# General configuration options:
#

OPTION(DEAL_II_WITH_FUNCTIONPARSER
  "Build deal.II with support for functionparser."
  ON)

OPTION(DEAL_II_WITH_METIS
  "Build deal.II with support for Metis."
  ON)

OPTION(DEAL_II_WITH_MPI
  "Build deal.II with support for mpi."
  OFF)

OPTION(DEAL_II_WITH_NETCDF
  "Build deal.II with support for netcdf."
  ON)

OPTION(DEAL_II_WITH_THREADS
  "Build deal.II with support for threads. This pulls in libtbb as a dependency."
  ON)

OPTION(DEAL_II_WITH_UMFPACK
  "Build deal.II with support for UMFPACK, BLAS and LAPACK."
  ON)

OPTION(DEAL_II_WITH_ZLIB
  "Build deal.II with support for zlib."
  ON)

#
# Options for library selection:
#

OPTION(DEAL_II_ALLOW_CONTRIB
  "Allow the use of contrib libraries bundled with the source tarball.
  Note: If set to off DEAL_II_FORCE_CONTRIB* will still pull in bundled
  packages, so to ensure that only external libraries are used
  DEAL_II_ALLOW_CONTRIB as well as all DEAL_II_FORCE_CONTRIB_* have to be
  OFF"
  ON)

OPTION(DEAL_II_FORCE_CONTRIB_FUNCTIONPARSER
  "Always use the bundled functionparser library instead of an external one."
  OFF)

OPTION(DEAL_II_FORCE_CONTRIB_BOOST
  "Always use the bundled boost library instead of an external one."
  OFF)

OPTION(DEAL_II_FORCE_CONTRIB_TBB
  "Always use the bundled tbb library instead of an external one."
  OFF)

OPTION(DEAL_II_FORCE_CONTRIB_UMFPACK
  "Always use the bundled umfpack library instead of an external one."
  OFF)

#
# Compatibility support:
#

OPTION(DEAL_II_COMPAT_MAPPING
  "Backward compatibility support for functions and classes that do not
  take an explicit mapping variable, but rather use a default Q1 mapping
  instead."
  ON)



#
# General configuration for cmake:
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

SET(CMAKE_INCLUDE_CURRENT_DIR true)

SET(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/contrib/cmake/check/"
  "${CMAKE_SOURCE_DIR}/contrib/cmake/configure/"
  "${CMAKE_SOURCE_DIR}/contrib/cmake/macros/"
  "${CMAKE_SOURCE_DIR}/contrib/cmake/modules/"
  )




#
# Set up deal.II specific variables
#

SET(DEAL_II_PACKAGE_BUGREPORT "dealii@dealii.org")
SET(DEAL_II_PACKAGE_NAME "deal.II")
SET(DEAL_II_PACKAGE_TARNAME ${DEAL_II_PACKAGE_NAME})
SET(DEAL_II_PACKAGE_VERSION ${VERSION})
SET(DEAL_II_PACKAGE_STRING "${DEAL_II_PACKAGE_NAME} ${DEAL_II_PACKAGE_VERSION}")
SET(DEAL_II_PACKAGE_URL "")
STRING(REGEX REPLACE "^([0-9]+)\\..*" "\\1" DEAL_II_MAJOR "${VERSION}")
STRING(REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1" DEAL_II_MINOR "${VERSION}")
SET(DEAL_II_PATH ${CMAKE_INSTALL_PREFIX})




#
# Now, configure a lot of things:
#

# used as a dummy target for all the stuff that has to be done before the
# library can be compiled:
ADD_CUSTOM_TARGET(deal_ii_target_dependencies)

# used to keep track of external libraries and necessary linker flags
SET(deal_ii_external_libraries "")
SET(deal_ii_external_debug_libraries "")
SET(deal_ii_required_linker_flags "")




#
# Feature configuration:
#

INCLUDE(configure_compiler_flags)

INCLUDE(configure_boost)

IF(DEAL_II_WITH_FUNCTIONPARSER)
  INCLUDE(configure_functionparser)
ENDIF()

IF(DEAL_II_WITH_MPI)
  INCLUDE(configure_mpi)
ENDIF()

IF(DEAL_II_WITH_NETCDF)
  INCLUDE(configure_netcdf)
ENDIF()

IF(DEAL_II_WITH_THREADS)
  INCLUDE(configure_threads)
ENDIF()

IF(DEAL_II_WITH_UMFPACK)
  INCLUDE(configure_umfpack)
ENDIF()

IF(DEAL_II_WITH_ZLIB)
  INCLUDE(configure_zlib)
ENDIF()




#
# Platform checks:
#

INCLUDE(check_for_compiler_features)

INCLUDE(check_for_cxx_features)

INCLUDE(check_for_compiler_bugs)




#
# Build the toolchain:
#

ADD_SUBDIRECTORY(common)

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/Version.in
  ${CMAKE_CURRENT_BINARY_DIR}/Version
  )




#
# Compile the library:
#

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_BINARY_DIR}/include/
  )

INCLUDE(macro_expand_instantiations)

ADD_SUBDIRECTORY(include)

ADD_SUBDIRECTORY(source)
