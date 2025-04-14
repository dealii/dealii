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
# A macro for the inst.in file expansion
#
# Usage:
#     expand_instantiations(target inst_in_files)
#
# Options:
#
# target
#
#    where target_${build_type} (and if present) target_${build_type}_cuda
#    will depend on the generation of all .inst files, to ensure that all
#    .inst files are generated prior to compiling.
#
# inst_in_files
#
#    a list of inst.in files that will be expanded
#

macro(expand_instantiations _target _inst_in_files)

  foreach (_inst_in_file ${_inst_in_files})
    string(REGEX REPLACE "\\.in$" "" _inst_file "${_inst_in_file}" )

    if(NOT CMAKE_CROSSCOMPILING)
      set(_command expand_instantiations_exe)
      set(_dependency expand_instantiations_exe)
    else()
      set(_command expand_instantiations)
      set(_dependency)
    endif()

    # create a .inst.tmp file first and only move to the correct name if the
    # first call succeeds. Otherwise we might be generating an incomplete
    # .inst file
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}
      DEPENDS ${_dependency}
              ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/template-arguments
              ${CMAKE_CURRENT_SOURCE_DIR}/${_inst_in_file}
      COMMAND ${_command}
      ARGS ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/template-arguments
           < ${CMAKE_CURRENT_SOURCE_DIR}/${_inst_in_file}
           > ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}.tmp
      COMMAND ${CMAKE_COMMAND}
      ARGS -E rename
           ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}.tmp
           ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}
      )

    list(APPEND _inst_targets ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file})
  endforeach()

  #
  # Define a custom target that depends on the generation of all inst.in
  # files.
  #
  add_custom_target(${_target}_inst ALL DEPENDS ${_inst_targets})

  #
  # Provide a way to generate all .inst files with a custom target.
  #
  add_dependencies(expand_all_instantiations ${_target}_inst)

  #
  # Add a dependency to all target.${build_type} so that target.inst is
  # fully generated before target will be processed.
  #
  foreach(_build ${DEAL_II_BUILD_TYPES})
    string(TOLOWER ${_build} _build_lowercase)

    add_dependencies(${_target}_${_build_lowercase} ${_target}_inst)

    if(TARGET ${_target}_${_build_lowercase}_cuda)
      add_dependencies(${_target}_${_build_lowercase}_cuda ${_target}_inst)
    endif()

  endforeach()

endmacro()
