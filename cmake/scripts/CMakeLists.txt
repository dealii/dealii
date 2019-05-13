## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2019 by the deal.II authors
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

ADD_EXECUTABLE(expand_instantiations_exe expand_instantiations.cc)

#
# We run into a bug when compiling expand_instantiations with
# -fuse-ld=gold -pthreads, see https://github.com/dealii/dealii/issues/1798
# Work around this by stripping -fuse-ld=gold out of the build flags
# for the script. This does little harm because linking this one file
# is so exceedingly cheap that there is no speed difference between
# using gold or the old BFD ld linker.
#
# If we are on a system where -fuse-ld=gold simply isn't part of the
# linker flag, filtering this one command out of the list of flags
# does not harm either.
#
STRING(REPLACE "-fuse-ld=gold" "" _expand_instantiations_link_flags "${DEAL_II_LINKER_FLAGS}" )

SET_TARGET_PROPERTIES(expand_instantiations_exe PROPERTIES
  LINK_FLAGS "${_expand_instantiations_link_flags}"
  LINKER_LANGUAGE "CXX"
  COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS}"
  COMPILE_FLAGS "${DEAL_II_CXX_FLAGS}"
  OUTPUT_NAME expand_instantiations
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}
  )

IF (DEAL_II_MSVC)
  SET_TARGET_PROPERTIES(expand_instantiations_exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
    )
ENDIF()

SET(_scripts
  ${CMAKE_CURRENT_SOURCE_DIR}/normalize.pl
  ${CMAKE_CURRENT_SOURCE_DIR}/run_test.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  ${CMAKE_CURRENT_SOURCE_DIR}/get_closest_tag.sh
  ${CMAKE_CURRENT_SOURCE_DIR}/get_latest_tag.sh
  )
FILE(COPY ${_scripts}
  DESTINATION ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/scripts
  )
INSTALL(FILES ${_scripts}
  DESTINATION ${DEAL_II_SHARE_RELDIR}/scripts
  COMPONENT library
  )
