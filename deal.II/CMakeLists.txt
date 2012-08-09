PROJECT(deal.II)
SET(VERSION "8.0.pre")
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

#TODO:
SET(CMAKE_CXX_FLAGS "-std=c++0x -D_REENTRANT -fPIC -O2 -march=native") #TODO: c++0x Features umschalten
SET(CMAKE_INCLUDE_CURRENT_DIR true)


#
# General configuration options:
#

OPTION(DEAL_II_USE_CONTRIB "Build and use contrib libraries bundled with the source tarball." OFF)
OPTION(DEAL_II_WITH_UMFPACK "Build deal.II with support for UMFPACK, BLAS and LAPACK." ON)
OPTION(DEAL_II_WITH_THREADS "Build deal.II with support for threads. This pulls in libtbb as a dependency." ON)
OPTION(DEAL_II_WITH_METIS "Build deal.II with support for Metis." OFF)
OPTION(DEAL_II_MPI "Build deal.II with support for openmpi." OFF)
OPTION(DEAL_II_FUNCTIONPARSER "BUild deal.II with support for functionparser" ON)
# TODO: The rest ;-)



#
# Configure:
#

SET(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/contrib/cmake/modules/"
  )

#
# DEAL_II_WITH_THREADS:
#

IF(DEAL_II_WITH_THREADS)
  FIND_PACKAGE(Threads REQUIRED)

  IF(DEAL_II_USE_CONTRIB)
    # Compiles and links libtbb, exports the TBB_* variables as well:
    ADD_SUBDIRECTORY(contrib/tbb)
  ELSE()
    FIND_PACKAGE(TBB REQUIRED)
  ENDIF()

  INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIR})

  SET(DEAL_II_USE_MT TRUE)
  SET(DEAL_II_USE_MT_POSIX TRUE)
ENDIF()


#
#
#


#
# Build the toolchain and contrib libraries:
#

ADD_SUBDIRECTORY(common/expand_instantiations)


#
# Compile the library:
#

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_BINARY_DIR}/include/
  )

INCLUDE(${CMAKE_SOURCE_DIR}/contrib/cmake/macros/macro_expand_instantiations.cmake)

ADD_SUBDIRECTORY(include)
ADD_SUBDIRECTORY(source)
