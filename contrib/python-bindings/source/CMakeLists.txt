## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2016 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

include_directories(
  ${CMAKE_BINARY_DIR}/include/
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_SOURCE_DIR}/contrib/python-bindings/include/
  SYSTEM
  ${DEAL_II_INCLUDE_DIRS}
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  )

#
# Build a variant of the wrapper library for all configured build types in
# DEAL_II_BUILD_TYPES (usually PyDealIIDebug and PyDealIIRelease)
#

include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

set(_src
  wrappers.cc
  export_tria_accessor.cc
  export_cell_accessor.cc
  export_point.cc
  export_triangulation.cc
  export_mapping.cc
  export_manifold.cc
  export_quadrature.cc
  cell_accessor_wrapper.cc
  tria_accessor_wrapper.cc
  point_wrapper.cc
  triangulation_wrapper.cc
  mapping_wrapper.cc
  manifold_wrapper.cc
  quadrature_wrapper.cc
  reference_cell_wrapper.cc
  )

foreach(_build ${DEAL_II_BUILD_TYPES})
  string(TOLOWER ${_build} _build_lowercase)

  PYTHON_ADD_MODULE(PyDealII_${_build_lowercase} ${_src})

  set_target_properties(PyDealII_${_build_lowercase} PROPERTIES
    OUTPUT_NAME "${PYTHON_BINDINGS_${_build}_NAME}"
    LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
    LINKER_LANGUAGE "CXX"
    COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
    COMPILE_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_PYTHON_RELDIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_PYTHON_RELDIR}"
    )
  if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set_target_properties(PyDealII_${_build_lowercase} PROPERTIES
      MACOSX_RPATH OFF
      BUILD_WITH_INSTALL_RPATH OFF
      INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${DEAL_II_PYTHON_RELDIR}"
      )
  endif()

  target_link_libraries(PyDealII_${_build_lowercase}
    ${DEAL_II_TARGET_NAME}_${_build_lowercase}
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
    )

  export(TARGETS PyDealII_${_build_lowercase}
    FILE ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
    APPEND
    )

  install(TARGETS  PyDealII_${_build_lowercase}
    COMPONENT python_bindings
    EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
    LIBRARY DESTINATION ${DEAL_II_PYTHON_RELDIR}
    ARCHIVE DESTINATION ${DEAL_II_PYTHON_RELDIR}
    )

  add_dependencies(python_bindings PyDealII_${_build_lowercase})
endforeach()

#
# Copy python sources to binary directory and set up installation:
#

set(_python_sources
  __init__.py
  )
file(COPY ${_python_sources}
  DESTINATION ${CMAKE_BINARY_DIR}/${DEAL_II_PYTHON_RELDIR}
  )
install(FILES ${_python_sources}
  DESTINATION ${DEAL_II_PYTHON_RELDIR}
  COMPONENT python_bindings
  )
