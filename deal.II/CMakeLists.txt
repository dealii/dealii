PROJECT(deal.II)

SET(VERSION "8.0.pre")

SET(DEAL_II_PACKAGE_BUGREPORT "dealii@dealii.org")
SET(DEAL_II_PACKAGE_NAME "deal.II")
SET(DEAL_II_PACKAGE_TARNAME ${DEAL_II_PACKAGE_NAME})
SET(DEAL_II_PACKAGE_VERSION ${VERSION})
SET(DEAL_II_PACKAGE_STRING "${DEAL_II_PACKAGE_NAME} ${DEAL_II_PACKAGE_VERSION}")
SET(DEAL_II_PACKAGE_URL "")




#
# General configuration for cmake:
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

SET(CMAKE_INCLUDE_CURRENT_DIR true)

SET(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/contrib/cmake/modules/"
  "${CMAKE_SOURCE_DIR}/contrib/cmake/macros/"
  "${CMAKE_SOURCE_DIR}/contrib/cmake/tests/"
  )




#
# General configuration options:
#

#TODO:
SET(CMAKE_CXX_FLAGS "-std=c++0x -Wfatal-errors -D_REENTRANT -fPIC -O2 -march=native") #TODO: c++0x Features umschalten

OPTION(DEAL_II_WITH_THREADS
  "Build deal.II with support for threads. This pulls in libtbb as a dependency."
  ON)

OPTION(DEAL_II_FUNCTIONPARSER
  "Build deal.II with support for functionparser."
  ON)

OPTION(DEAL_II_MPI
  "Build deal.II with support for openmpi."
  OFF)

OPTION(DEAL_II_WITH_UMFPACK
  "Build deal.II with support for UMFPACK, BLAS and LAPACK."
  ON)

OPTION(DEAL_II_WITH_METIS
  "Build deal.II with support for Metis."
  OFF)

OPTION(DEAL_II_ALLOW_CONTRIB
  "Allow the use of contrib libraries bundled with the source tarball.
  Note: If set to off DEAL_II_FORCE_CONTRIB* will still pull in bundled
  packages, so to ensure that only external libraries are used
  DEAL_II_ALLOW_CONTRIB as well as all DEAL_II_FORCE_CONTRIB_* have to be
  OFF"
  ON)

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
# Now, configure and configure a lot of things:
#

# used as a dummy target for all the stuff that has to be done before the
# library can be compiled:
ADD_CUSTOM_TARGET(deal_ii_target_dependencies)




#
# Platform checks:
#

INCLUDE(check_for_cxx_features)




#
# Feature configuration:
#

INCLUDE(configure_boost)

IF(DEAL_II_WITH_THREADS)
  INCLUDE(configure_threads)
ENDIF()

IF(DEAL_II_WITH_UMFPACK)
  INCLUDE(configure_umfpack)
ENDIF()




#
# Build the toolchain:
#

ADD_SUBDIRECTORY(common/expand_instantiations)




#
# Compile the library:
#

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_SOURCE_DIR}/contrib/ #TODO
  ${CMAKE_BINARY_DIR}/include/
  )

INCLUDE(macro_expand_instantiations)

ADD_SUBDIRECTORY(include)

ADD_SUBDIRECTORY(source)
