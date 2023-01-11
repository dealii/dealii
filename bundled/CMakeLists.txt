## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Compile and install enabled bundled features:
#

message(STATUS "Setting up bundled features")

include_directories(
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  ${DEAL_II_INCLUDE_DIRS}
  )

if(DEAL_II_FEATURE_BOOST_BUNDLED_CONFIGURED)
  add_subdirectory(${BOOST_FOLDER}/libs/serialization/src)
  add_subdirectory(${BOOST_FOLDER}/libs/system/src)

  install(DIRECTORY ${BOOST_FOLDER}/include/boost
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    )

  if(DEAL_II_WITH_ZLIB)
    add_subdirectory(${BOOST_FOLDER}/libs/iostreams/src)
  else()
    message(STATUS "BOOST::Iostreams will not be available because zlib or its header files could not be found")
  endif()
endif()


if(DEAL_II_FEATURE_KOKKOS_BUNDLED_CONFIGURED)
  add_subdirectory(${KOKKOS_FOLDER})

  install(DIRECTORY
    ${KOKKOS_FOLDER}/algorithms/src/
    ${KOKKOS_FOLDER}/containers/src/
    ${KOKKOS_FOLDER}/core/src/
    ${KOKKOS_FOLDER}/simd/src/
    ${KOKKOS_FOLDER}/tpls/desul/include/
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    FILES_MATCHING PATTERN "*.hpp"
    )
endif()


if(DEAL_II_FEATURE_TASKFLOW_BUNDLED_CONFIGURED)
  install(DIRECTORY ${TASKFLOW_FOLDER}/include/taskflow
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    )
endif()


if(DEAL_II_FEATURE_TBB_BUNDLED_CONFIGURED)
  add_subdirectory(${TBB_FOLDER}/src)

  install(DIRECTORY ${TBB_FOLDER}/include/tbb
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    FILES_MATCHING PATTERN "*.h"
    )
endif()


if(DEAL_II_FEATURE_UMFPACK_BUNDLED_CONFIGURED)
  add_subdirectory(${UMFPACK_FOLDER}/UMFPACK/Source)
  add_subdirectory(${UMFPACK_FOLDER}/AMD/Source)

  install(DIRECTORY
      ${UMFPACK_FOLDER}/UMFPACK/Include/
      ${UMFPACK_FOLDER}/AMD/Include/
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    FILES_MATCHING PATTERN "*.h"
    )
endif()


if(DEAL_II_FEATURE_MUPARSER_BUNDLED_CONFIGURED)
  add_subdirectory(${MUPARSER_FOLDER})
endif()


message(STATUS "Setting up bundled features - Done")
