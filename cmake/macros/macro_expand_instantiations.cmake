## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# A macro for the inst.in file expansion
#
# Usage:
#     EXPAND_INSTANTATIONS(target inst_in_files)
#
# Options:
#
# target
#
#    where target.${build_type} will depend on the generation of all .inst
#    files, to ensure that all .inst files are generated prior to
#    compiling.
#
# inst_in_files
#
#    a list of inst.in files that will be expanded
#

MACRO(EXPAND_INSTANTIATIONS _target _inst_in_files)

  FOREACH (_inst_in_file ${_inst_in_files})
    STRING(REGEX REPLACE "\\.in$" "" _inst_file "${_inst_in_file}" )

    IF(NOT CMAKE_CROSSCOMPILING)
      SET(_command expand_instantiations_exe)
      SET(_dependency expand_instantiations_exe)
    ELSE()
      SET(_command expand_instantiations)
      SET(_dependency)
    ENDIF()

    # create a .inst.tmp file first and only move to the correct name if the
    # first call succeeds. Otherwise we might be generating an incomplete
    # .inst file
    ADD_CUSTOM_COMMAND(
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

    LIST(APPEND _inst_targets ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file})
  ENDFOREACH()

  #
  # Define a custom target that depends on the generation of all inst.in
  # files.
  #
  ADD_CUSTOM_TARGET(${_target}.inst ALL DEPENDS ${_inst_targets})

  #
  # Add a dependency to all target.${build_type} so that target.inst is
  # fully generated before target will be processed.
  #
  FOREACH(_build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${_build} _build_lowercase)
    ADD_DEPENDENCIES(${_target}.${_build_lowercase} ${_target}.inst)
  ENDFOREACH()

ENDMACRO()

