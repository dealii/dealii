## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# Compile and install enabled bundled features:
#

message(STATUS "Setting up bundled features")

#
# If we install any bundled software also install the README.md file into
# the bundled header directory. First of all, this ensures that we always
# have a disclaimer in place explaining where the header files come from.
# Secondly, in rare circumstances we might end up not installing any header
# (because they are only used internally). Because this is quite tricky to
# figure out when setting up targets we work around this issue by simply
# always installing something into the bundled header directory.
#
if(NOT "${DEAL_II_BUNDLED_INCLUDE_DIRS}" STREQUAL "")
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/README.md
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    )
endif()

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
    ${KOKKOS_FOLDER}/tpls/mdspan/include/
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    FILES_MATCHING
    PATTERN "*.hpp"
    PATTERN "*.h"
    )
endif()


if(DEAL_II_FEATURE_TASKFLOW_BUNDLED_CONFIGURED)
  install(DIRECTORY ${TASKFLOW_FOLDER}/taskflow
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
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
