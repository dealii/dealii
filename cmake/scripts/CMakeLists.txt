## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

add_executable(expand_instantiations_exe expand_instantiations.cc)

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
string(REPLACE "-fuse-ld=gold" "" _expand_instantiations_link_flags "${DEAL_II_LINKER_FLAGS}" )

set_target_properties(expand_instantiations_exe PROPERTIES
  LINK_FLAGS "${_expand_instantiations_link_flags}"
  LINKER_LANGUAGE "CXX"
  COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS}"
  COMPILE_FLAGS "${DEAL_II_CXX_FLAGS}"
  OUTPUT_NAME expand_instantiations
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}
  )

if (DEAL_II_MSVC)
  set_target_properties(expand_instantiations_exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
    )
endif()

set(_scripts
  ${CMAKE_CURRENT_SOURCE_DIR}/normalize.pl
  ${CMAKE_CURRENT_SOURCE_DIR}/run_test.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  ${CMAKE_CURRENT_SOURCE_DIR}/get_closest_tag.sh
  ${CMAKE_CURRENT_SOURCE_DIR}/get_latest_tag.sh
  ${CMAKE_CURRENT_SOURCE_DIR}/get_fancy_tag.sh
  )
file(COPY ${_scripts}
  DESTINATION ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/scripts
  )
install(FILES ${_scripts}
  DESTINATION ${DEAL_II_SHARE_RELDIR}/scripts
  COMPONENT library
  )
