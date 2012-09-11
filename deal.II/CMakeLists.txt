PROJECT(deal.II)

SET(VERSION "8.0.pre")




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



#TODO:
SET(CMAKE_CXX_FLAGS "-std=c++0x -Wfatal-errors -D_REENTRANT -fPIC -O2 -march=native")
SET(CMAKE_C_FLAGS "-Wfatal-errors -D_REENTRANT -fPIC -O2 -march=native")




#
# General configuration options:
#

OPTION(DEAL_II_WITH_FUNCTIONPARSER
  "Build deal.II with support for functionparser."
  ON)

OPTION(DEAL_II_WITH_METIS
  "Build deal.II with support for Metis."
  OFF)

OPTION(DEAL_II_WITH_MPI
  "Build deal.II with support for mpi."
  ON)

OPTION(DEAL_II_WITH_NETCDF
  "Build deal.II with support for netcdf."
  ON)

OPTION(DEAL_II_WITH_THREADS
  "Build deal.II with support for threads. This pulls in libtbb as a dependency."
  ON)

OPTION(DEAL_II_WITH_UMFPACK
  "Build deal.II with support for UMFPACK, BLAS and LAPACK."
  ON)



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
# Now, configure a lot of things:
#

# used as a dummy target for all the stuff that has to be done before the
# library can be compiled:
ADD_CUSTOM_TARGET(deal_ii_target_dependencies)

# used to keep track of external libraries and necessary link flags
SET(deal_ii_external_libraries "")
SET(deal_ii_external_debug_libraries "")
SET(deal_ii_required_link_flags "")




#
# Feature configuration:
#

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





#
# Platform checks:
#

INCLUDE(check_for_cxx_features)

INCLUDE(check_for_compiler_bugs)




#
# Build the toolchain:
#

ADD_SUBDIRECTORY(common/expand_instantiations)




#
# Compile the library:
#

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/include/
  # ${CMAKE_SOURCE_DIR}/contrib/ #TODO
  ${CMAKE_BINARY_DIR}/include/
  )

INCLUDE(macro_expand_instantiations)

ADD_SUBDIRECTORY(include)

ADD_SUBDIRECTORY(source)

