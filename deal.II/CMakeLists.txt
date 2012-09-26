
SET(VERSION "8.0.pre")

###########################################################################
##                                                                       ##
##  See README.cmake for further details on how to use the cmake build   ##
##  system of deal.II.                                                   ##
##                                                                       ##
###########################################################################

#
#             The cmake build system for the deal.II project
#
# The main purpose of this file is to set all the configuration options and
# then call a lot of modules residing under ./cmake/*/ that do the
# hard work.
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
# Enables feature autodetection. This will automatically overwrite all
# DEAL_II_WITH_<...> toggles depending on whether they can be supported or
# not.
# Short description:
  "Enables feature autodetection. This will automatically set DEAL_II_WITH_<...> options."
  ON)

OPTION(DEAL_II_ALLOW_CONTRIB
#  Allow the use of contrib libraries bundled with the source tarball.
#  Note: If set to off DEAL_II_FORCE_CONTRIB* will still pull in bundled
#  packages, so to ensure that only external libraries are used
#  DEAL_II_ALLOW_CONTRIB as well as all DEAL_II_FORCE_CONTRIB_* have to be
#  OFF
# Short description:
  "Allow the use of contrib libraries bundled with the source tarball. (DEAL_II_FORCE_CONTRIB* will overwrite this option.)"
  OFF)

OPTION(DEAL_II_WITH_ARPACK
  "Build deal.II with support for ARPACK."
  OFF)

OPTION(DEAL_II_WITH_BLAS
  "Build deal.II with support for BLAS."
  OFF)

SET(DEAL_II_WITH_BOOST
  # Always true. We need it :-]
  ON
  CACHE BOOL "Build deal.II with support for boost." FORCE
  )
OPTION(DEAL_II_FORCE_CONTRIB_BOOST
  "Always use the bundled boost library instead of an external one."
  OFF)

OPTION(DEAL_II_WITH_DOXYGEN
  "Build deal.II with support for doxygen and dot."
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
  "Build deal.II with support for p4est."
  OFF)

OPTION(DEAL_II_WITH_PETSC
  "Build deal.II with support for petsc."
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
#OPTION(DEAL_II_FORCE_CONTRIB_UMFPACK
#  "Always use the bundled umfpack library instead of an external one."
#  OFF)

OPTION(DEAL_II_WITH_ZLIB
  "Build deal.II with support for zlib."
  OFF)

#
# Build configuration: configuration options regarding compilation and
# installation of the deal.II library
#

OPTION(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS
  "configure sensible default CFLAGS and CXXFLAGS depending on platform, compiler and build target."
  ON)

#
# Component selection: configuration options regarding the setup of
# different components of the deal.II library
#

OPTION(DEAL_II_COMPONENT_EXAMPLES
  "Enable configuration and installation of the example steps. This adds a COMPONENT \"examples\" to the build system."
  ON)

OPTION(DEAL_II_COMPONENT_DOCUMENTATION
  "Enable configuration, build and installation of the documentation. This adds a COMPONENT \"documentation\" to the build system."
  OFF)

OPTION(DEAL_II_COMPONENT_COMPAT_FILES
  "Enable installation of the example steps. This adds a COMPONENT \"compat_files\" to the build system."
  ON)

OPTION(DEAL_II_COMPONENT_PROJECT_CONFIG
  # Enable configuration and installation of a cmake project config.
  # Useful for autodetection ( FIND_PACKAGE(DEAL_II) ) of the deal.II
  # library out of cmake build systems.
  "Enable configuration and installation of a cmake project config. This adds a COMPONENT \"project_config\" to the build system."
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

LIST(APPEND CMAKE_MODULE_PATH
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


#
# Run all system checks:
#
FILE(GLOB check_files "cmake/checks/*.cmake")
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

