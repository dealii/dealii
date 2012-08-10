PROJECT(deal.II)
SET(VERSION "8.0.pre")
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

#TODO:
SET(CMAKE_CXX_FLAGS "-std=c++0x -D_REENTRANT -fPIC -O2 -march=native") #TODO: c++0x Features umschalten
SET(CMAKE_INCLUDE_CURRENT_DIR true)



#
# General configuration options:
#

OPTION(DEAL_II_USE_CONTRIB
  "Build and use contrib libraries bundled with the source tarball."
  ON
  )
OPTION(DEAL_II_WITH_THREADS
  "Build deal.II with support for threads. This pulls in libtbb as a dependency."
  ON
  )
OPTION(DEAL_II_FUNCTIONPARSER
  "BUild deal.II with support for functionparser"
  ON
  )
OPTION(DEAL_II_MPI "Build deal.II with support for openmpi."
  OFF
  )
OPTION(DEAL_II_WITH_UMFPACK "Build deal.II with support for UMFPACK, BLAS and LAPACK."
  ON
  )
OPTION(DEAL_II_WITH_METIS "Build deal.II with support for Metis."
  OFF
  )

# TODO: The rest ;-)



#
# Now, configure and prepare a lot of things:
#

SET(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/contrib/cmake/modules/"
  )

#
# used as a dummy target for all the stuff that has to be done before the
# library can be compiled:
#
ADD_CUSTOM_TARGET(deal_ii_target_dependencies)


# TODO: Version Mumbo Jumbo

IF(DEAL_II_WITH_THREADS)
  INCLUDE(preprare_threads)
ENDIF()

IF(DEAL_II_WITH_UMFPACK)
  INCLUDE(prepare_umfpack)
ENDIF()

# TODO: The rest...

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
